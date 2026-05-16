# SymPy 実装ガイドライン — dppu エンジン

- **Target:** AI Coding Assistants & Contributors
- **Context:** Einstein-Cartan Theory on Curved Spacetime (S³×S¹, T³×S¹, Nil³×S¹)

本ドキュメントは、`dppu` パッケージにおける数式処理（SymPy）の実装指針、および理論物理学的規約（Convention）を定めるものである。以下のルールを厳守すること。

⇒ [English version](SymPy_guideline.md) | [幾何規約](CONVENTIONS_ja.md)

-----

## 1\. 高速化エンジニアリング指針 (Optimization Rules)

SymPyを用いた積分計算における「処理時間爆発（2時間以上）」を防ぎ、数秒で完了させるための鉄則。

### Rule 1.1: `expand()` + `cancel()` 戦略の徹底

積分（`integrate`）を行う直前に、被積分関数に対して高コストな `simplify()` を使用してはならない。代わりに「展開と約分」を行うこと。

  * **Don't:**
    ```python
    density = simplify(density)  # NG 禁止：巨大な式の因数分解は計算コストが極大
    result = integrate(density, x)
    ```
  * **Do:**
    ```python
    density = cancel(expand(density))  # OK 推奨：多項式の和にすることで項別積分を誘発
    result = integrate(density, x)
    ```

### Rule 1.2: 積分後の中間簡約の抑制

多重積分（例： $\phi$ 積分 → $\theta$ 積分）の際、中間結果に対して過度な `simplify` を行わない。`cancel` 程度に留め、最終的な積分直前で再度 `expand` する方が高速である。

### Rule 1.3: Mixing DOF のための cos/sin 多項式補助シンボル

mixing 回転角 $\theta_i = \arctan\delta_i$ が構造定数に現れる場合（fiber mode MIXING または BOTH）、パイプライン中に `1/sqrt(1+delta_i**2)` を直接書いてはならない。代わりに**補助シンボル** `cd_i`、`sd_i` を使用する：

$$
cd_i \equiv \cos\theta_i = \frac{1}{\sqrt{1+\delta_i^2}},\quad
sd_i \equiv \sin\theta_i = \frac{\delta_i}{\sqrt{1+\delta_i^2}}.
$$

エンジンは `_build_fiber_params`（MIXING/BOTH モード時）でこれらを自動的に導入する。全ての中間式が `{r, L, eta, V, omega_k, cd_i, sd_i, ...}` の**多項式**として扱われ、E4 パイプライン全体で `expand`/`cancel` が高速に機能する。

**なぜ必要か**：繰り返し現れる `1/sqrt(1+δ²)` は SymPy に多項式構造を認識させない。`cancel` が計算上扱いにくい形を生成し、実行時間が時間単位になる。`cd_i` 置換で分単位に短縮できる。

**数値代入**（ $\delta=0$、等方 mixing 点）：
```python
iso_subs = {fp['cd0']: 1, fp['sd0']: 0,
            fp['cd1']: 1, fp['sd1']: 0,
            fp['cd2']: 1, fp['sd2']: 0}
Veff_numeric = Veff_sym.subs(iso_subs).subs(phys_subs)
```

**$\delta_i$ に関する微分**（ $\delta=0$）： $\partial/\partial\delta_i\big|\_0 = \partial/\partial(sd_i)\big|\_{sd=0,\,cd=1}$ を利用し、`sd_i` で微分後 `sd_i=0, cd_i=1` を代入する：
```python
dV_ddelta0 = diff(Veff_sym, fp['sd0']).subs(iso_subs)
```

### Rule 1.4: 数値評価のための `lambdify` + CSE

数値 Hessian 計算やパラメータスキャンには、シンボル式 $V_{\rm eff}$ を `lambdify` で高速 NumPy 関数に変換する。CSE（共通部分式除去）を有効化することで、数千項の式でも大幅に高速化できる：

```python
from sympy import lambdify

fp = engine.get_free_params()          # 有効な Symbol の dict
sym_args = tuple(fp[k] for k in ordered_keys)

Vfunc = lambdify(sym_args, Veff_expr, modules='numpy', cse=True)
# cse=True: 数千項の式で大幅高速化
```

S³ BOTH モードの典型的な `ordered_keys`：
```python
ordered_keys = ['r', 'L', 'kappa', 'eta', 'V', 'theta_NY', 'alpha',
                'omega1', 'omega2', 'omega3',
                'cd0', 'sd0', 'cd1', 'sd1', 'cd2', 'sd2']
```

等方評価点では `cd_i=1, sd_i=0`（全 $i$）を渡す。

### Rule 1.5: 大きなシンボル行列に対して `Matrix.inv()` を使わない

複数のシンボルを含む SymPy `Matrix` に `.inv()` を呼ぶと、無期限のハングが発生することがある（ $4\times4$ 多項式成分行列で >30 分の観測例あり）。

* **NG**：`G_inv = G_matrix.inv()`
* **OK**：小行列には明示的な閉形式逆行列を用いるか、`lambdify` 後に数値的に逆行列を計算する。

```python
# 例：明示的 2×2 逆行列
det = a * d - b * c
G_inv = Matrix([[d, -b], [-c, a]]) / det
```

この問題は mixing フレームの構造定数変換で発生した。`dppu/topology/s3.py` で明示的 `G_inv` 公式に修正済み。

-----

## 2\. 理論実装指針 (Theoretical Implementation Rules)

曲がった時空（Curved Spacetime）において、テンソル演算の整合性を保つための鉄則。

### Rule 2.1: 配列インデックス操作の禁止 (Robust Method)

非対角計量や $g_{\mu\nu} \neq 1$ の環境下では、`T[mu, nu, lam]` のように配列のインデックスを入れ替える操作は、物理的なテンソルの添字操作（上げ下げ）と等価ではない。

  * **Don't:**
    ```python
    # NG 禁止：物理的に誤った成分を参照する恐れがある
    term = T_tensor[mu, nu, lam]
    ```
  * **Do:** 必ず計量 $g_{\mu\nu}$ を介して操作する。
    1.  全ての添字を下げて完全共変形 $T_{\lambda\mu\nu}$ を作る。
    2.  添字の入れ替え（Permutation）を行う。
    3.  必要に応じて計量で添字を上げる。

#### 正規直交フレーム基底での最適化

正規直交フレーム基底においては、計量が単位行列（ $g_{ab} = \eta_{ab} = \text{diag}(1,1,1,1)$ または $\text{diag}(-1,1,1,1)$）であるため、添字の上げ下げ計算を省略し、直接成分計算を行うことで高速化する。

ただし、物理的定義（Hehl 1976）の符号パターン $(+1, +1, -1)$ は厳守すること。

Lorentzian/LZ-native 経路ではこの省略を行わない。計量は $\mathrm{diag}(-1,+1,+1,+1)$ であり、`contortion`、Weyl scalar、Pontryagin/Hodge contraction は必ず `metric` / `metric_inv` による添字操作を使う。下のコード断片は legacy Euclidean 経路または恒等計量の説明用である。

```python
# ============================================================
# Golden Logic for Contortion (Frame Basis / legacy Euclidean path)
# ============================================================
# Assumption: Metric is diagonal/identity (Orthonormal Frame)
# Therefore T^a_bc and T_abc behave identically in code logic.

K_tensor = MutableDenseNDimArray.zeros(dim, dim, dim)

for a in range(dim):
    for b in range(dim):
        for c in range(dim):
            # Formula: K_abc = (1/2)(T_abc + T_bca - T_cab)
            # Note: Using T[a,b,c] directly as T_abc
            
            term = (T_tensor[a, b, c] + T_tensor[b, c, a] - T_tensor[c, a, b])
            
            val = term * Rational(1, 2)
            
            if val != 0:
                K_tensor[a, b, c] = cancel(expand(val))
```

### Rule 2.2: 自己無撞着性の検証 (Consistency Check)

EC接続 ($\Gamma_{\text{EC}}$) を構築した後、必ず以下の検証コードを実行し、ミスマッチが **0** であることを確認しなければならない。

```python
# Torsion Consistency Check
T_verify = Gamma_EC[lam, mu, nu] - Gamma_EC[lam, nu, mu] # Hehl定義
mismatch = count(simplify(T_verify - T_original) != 0)
assert mismatch == 0
```

-----

## 3\. 理論的規約 (Standard Conventions)

論文執筆時の混乱を防ぐため、**Hehl (1976) 標準** に準拠する。

### 3.1 Torsion Definition

捩れテンソル $T^\lambda_{\ \mu\nu}$ の定義：
$$T^\lambda_{\ \mu\nu} \equiv \Gamma^\lambda_{\ \mu\nu} - \Gamma^\lambda_{\ \nu\mu}$$
（注：本エンジンのフレーム基底での捩れ成分は、CONVENTIONS の第6節で定義された捩れ2-form $T^a = de^a + \omega^{a}{}\_b\wedge e^b$ から $T^a = \frac{1}{2}T^{a}{}\_{bc}\,e^b\wedge e^c$ の係数比較により抽出される）

### 3.2 Contortion Formula

上記のTorsion定義と整合するContortion $K^\lambda_{\ \mu\nu}$ の公式（Verified Formula）：

$$K_{\lambda\mu\nu} = \frac{1}{2} \left( T_{\lambda\mu\nu} + T_{\mu\nu\lambda} - T_{\nu\lambda\mu} \right)$$

  * 符号パターン: **$(+1, +1, -1)$**
  * 注意: この式は $T_{\lambda\mu\nu}$（全ての添字を下げたもの）に対して適用すること。

### 3.3 Einstein-Cartan Connection

$${\Gamma_{\text{EC}}}^\lambda_{\ \mu\nu} = {\Gamma_{\text{LC}}}^\lambda_{\ \mu\nu} + K^\lambda_{\ \mu\nu}$$

  * $\Gamma_{\text{LC}}$: Levi-Civita接続（Christoffel記号）
  * $K$: Contortion

-----

## 4. Torsion Ansatz と Mode 分解規約

$M^3 \times S^1$ ミニスーパースペース・アンサーツでの捩れテンソルは、以下の3モードで指定する。

### 4.1 モード定義

| Mode | 物理的成分 | パラメータ |
|---|---|---|
| `Mode.AX` | 軸対称成分（T1）のみ | $\eta \neq 0$, $V = 0$ |
| `Mode.VT` | ベクトル跡部分（T2）のみ | $\eta = 0$, $V \neq 0$ |
| `Mode.MX` | T1 + T2 の両成分 | $\eta \neq 0$, $V \neq 0$ |

### 4.2 物理的対応

- **T1（軸対称成分）**: 軸性ベクトル $S^\mu = (\eta/r)(0,0,0,1)$ の双対。空間添字 $a,b,c \in \{0,1,2\}$ に対して $T_{abc} = (2\eta/r)\,\varepsilon_{abc}$。
- **T2（ベクトル跡部分）**: ベクトル $V_\mu = V\,\delta^3_\mu$（$\tau$ 成分のみ）の双対。 $T_{abc} = \frac{1}{3}(\delta_{ac}V_b - \delta_{ab}V_c)$。

### 4.3 実装ルール

`dppu/torsion/ansatz.py` の `construct_torsion_tensor(mode, r, eta, V, metric, dim)` を使用すること。
$T_{abc}$ を手打ちで構築してはならない。

-----

## 5. Nieh-Yan トポロジカル項のバリアント

### 5.1 Nieh-Yan 分解

完全な Nieh-Yan 密度：
$$N = N_{\mathrm{TT}} - N_{\mathrm{Ree}},$$
$$N_{\mathrm{TT}} = \frac{1}{4}\varepsilon^{abcd}T^{e}{}\_{ab}T_{ecd},\qquad N_{\mathrm{Ree}} = \frac{1}{4}\varepsilon^{abcd}R_{abcd}.$$

### 5.2 バリアント選択

| `NyVariant` | 使用する密度 |
|---|---|
| `NyVariant.TT` | $N_{\mathrm{TT}}$ のみ |
| `NyVariant.REE` | $N_{\mathrm{Ree}}$ のみ |
| `NyVariant.FULL` | $N_{\mathrm{TT}} - N_{\mathrm{Ree}}$（標準） |

### 5.3 実装

パイプラインステップ `E4.10` で全バリアントが計算される。バリアント選択は engine の `__init__` で `ny_variant` 引数を指定する。`dppu/torsion/nieh_yan.py` を参照。

-----

## 6. 拡張ラグランジアンと Weyl 結合定数 $\alpha$

### 6.1 作用の形式

$$S = \int L \times \mathrm{Vol},\qquad
L = \frac{R}{2\kappa^2} + \theta_{\mathrm{NY}}\times N + \alpha\times C^2.$$

| パラメータ | 意味 |
|---|---|
| $\kappa$ | Einstein-Cartan 重力結合定数 |
| $\theta_{\mathrm{NY}}$ | Nieh-Yan 結合定数（トポロジカル） |
| $\alpha$ | Weyl 結合定数（共形不変項） |

$\alpha \leq 0$ では定理1により安定真空が保護される。 $\alpha > 0$ では定理2により有効ポテンシャルが非有界となる。

### 6.2 有効ポテンシャルの取得

**UnifiedEngine** ：
```python
engine.run()
fp   = engine.get_free_params()          # 有効な SymPy Symbol の dict
Veff = engine.data['potential']          # シンボル式

# 高速数値評価（Rule 1.4 参照）：
from sympy import lambdify
Vfunc = lambdify(tuple(fp[k] for k in ordered_keys), Veff, modules='numpy', cse=True)
```

$V_{\rm eff} = -S$。パイプラインステップ `E4.13` で抽出される。

### 6.3 実装位置

`dppu/action/lagrangian.py` の `compute_lagrangian()`、`dppu/action/potential.py` を参照。

-----

## 7. 数値最適化戦略（Phase Atlas 探索）

### 7.1 2段階戦略

**Stage 1：Brute-force グリッド探索**

`scipy.optimize.brute`（`Ns` 点/軸）で $(r, \varepsilon)$ 2D グリッドを粗く探索し、大域極小の領域を特定する。

**Stage 2：Multi-start L-BFGS-B 精密化**

Stage 1 の上位 $N$ 候補を初期点として `scipy.optimize.minimize`（L-BFGS-B、`ftol=1e-8`）を実行し、高精度最小値を得る。

### 7.2 安定性分類

$(r^\*, \varepsilon^\*)$ の最小値発見後、以下に分類する：

| 分類 | 条件 | 物理的意味 |
|---|---|---|
| Type-I | $V(r)$ が $r \to 0$ で立ち上がる（障壁あり） | 安定真空（核生成障壁あり） |
| Type-II | $V(r)$ が $r \to 0$ で単調減少（障壁なし） | 自発的核生成が可能 |
| Type-III | 物理的領域に局所最小なし | 不安定配位 |

`dppu/scanning/stability.py` の `analyze_stability()` を使用すること。

### 7.3 注意事項

- $\alpha > 0$ の場合、最適化が探索境界 $(r \to 0^+,\,\varepsilon \to -1^+)$ に張り付く（`converged = False`）。これは最適化の失敗ではなく、ポテンシャルが探索範囲内で非有界であることの正確な報告である。

-----

## 8. Self-duality（SD）診断規約

### 8.1 曲率の Hodge 双対

$$(*R)^{ab}{}_{cd} = \frac{1}{2}\varepsilon_{cdef}\,R^{ab,ef},$$

ここで $\varepsilon_{cdef}$ はフレーム基底の Levi-Civita 記号（`dppu/utils/epsilon.py` の `epsilon_4d()` を使用）。

### 8.2 Pontryagin 内積と SD 残差

$$E_{RR} = \langle R, R\rangle = R_{abcd}R^{abcd},\qquad
P = \langle R, *R\rangle = R_{abcd}(*R)^{abcd}.$$

| 条件 | 物理的状態 |
|---|---|
| SD 残差 $< \varepsilon_{\rm SD}$ かつ $\|R\| > \varepsilon_R$ | Self-dual |
| ASD 残差 $< \varepsilon_{\rm SD}$ かつ $\|R\| > \varepsilon_R$ | Anti-self-dual |
| `P_form = 0` | AX/VT/MX の全 EC branch で form-Hodge Pontryagin density は厳密にゼロ。 |
| `P_int = 0` — AX/VT | internal-pair Pontryagin diagnostic は AX/VT でゼロ |
| `P_int \neq 0` — MX | $P_{\rm int}=2V\eta(-V^2r^2+9\eta^2-36)/(9r^3)$。源は捩れ混合 $V\eta$ であり、twist ではない |

### 8.3 使用方法

```python
from dppu.curvature.sd_extension import SDExtensionMixin
SDExtensionMixin.attach_to(engine)          # メソッドを動的にアタッチ
R = engine.get_R_ab_cd_numerical(params_dict)
diag = engine.evaluate_sd_status(params_dict)
# diag['P_RstarR'] で数値 Pontryagin 診断値を確認
```

`dppu/curvature/sd_extension.py`、`dppu/curvature/sd_diagnostics.py`、`dppu/curvature/pontryagin.py` を参照。

-----

## 9\. 数値 Hessian 解析

### 9.1 基本パターン

等方点における有限差分 Hessian により安定性を評価する（記号的微分ではなく）。

```python
import numpy as np
from sympy import lambdify

# 1. エンジン実行・高速数値関数の構築（Rule 1.4）
fp = engine.get_free_params()
Veff_sym = engine.data['potential']
Vfunc = lambdify(tuple(fp[k] for k in ordered_keys), Veff_sym, modules='numpy', cse=True)

# 2. 摂動点で V_eff を評価するラッパー
phys_base = dict(r=r0, L=1, kappa=1, eta=1, V=1, theta_NY=0, alpha=0,
                 omega1=0, omega2=0, omega3=0,
                 cd0=1, sd0=0, cd1=1, sd1=0, cd2=1, sd2=0)

def veff(**overrides):
    vals = {**phys_base, **overrides}
    return float(Vfunc(*[vals[k] for k in ordered_keys]))

# 3. 中心差分 Hessian
h = 1e-4
keys = ['omega1', 'omega2', 'omega3', 'sd0', 'sd1', 'sd2']  # δ_i の微分は sd_i で行う
H = np.zeros((6, 6))
for i, ki in enumerate(keys):
    for j, kj in enumerate(keys):
        H[i, j] = (veff(**{ki: +h, kj: +h}) - veff(**{ki: +h, kj: -h})
                 - veff(**{ki: -h, kj: +h}) + veff(**{ki: -h, kj: -h})) / (4 * h**2)
```

### 9.2 $\delta_i$ に関する微分の規約

`delta_i` は `cd_i`, `sd_i` のみを通じて現れるため、有限差分変数は `sd_i`（`delta_i` ではない）を使う：

$$
\frac{\partial V}{\partial \delta_i}\bigg|_{\delta=0}
= \frac{\partial V}{\partial sd_i}\bigg|_{sd=0,\,cd=1}.
$$

有限差分呼び出しでは `sd_i = ±h`（`cd_i = 1` 固定）を渡す。

### 9.3 Spin-1 の一般化固有値問題

Spin-1 安定性の評価には $H v = \lambda G v$ を解く（`scipy.linalg.eigh`）：

```python
from scipy.linalg import eigh

# G は enable_velocity=True で計算（CONVENTIONS §12.4 参照）
# G は 6D で rank-3 → 先に物理部分空間に射影する

# null 方向 Y_k ∝ (1, 2L/r0)（各 k ∈ {0,1,2}）
norm_Y = np.sqrt(1 + (2*L/r0)**2)
Y_blocks = [np.array([1, 0, 0, 2*L/r0, 0, 0]) / norm_Y,   # k=0
            np.array([0, 1, 0, 0, 2*L/r0, 0]) / norm_Y,   # k=1
            np.array([0, 0, 1, 0, 0, 2*L/r0]) / norm_Y]   # k=2

# H, G を X_k（Y_k の直交補空間）に射影して 3×3 行列を得る

eigenvalues, _ = eigh(H_phys, G_phys)
# lambda_phys > 0  ⟹  物理 spin-1 セクターにタキオンなし
```

> **警告**：`eigh(H_full, G_full)` を直接呼ぶと `G` が特異のため `scipy` が `LinAlgError` を発生させる。必ず先に物理部分空間に射影すること。
