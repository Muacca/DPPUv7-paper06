# CONVENTIONS — dppu エンジン幾何規約

本書は、`BaseFrameEngine` と各 topology runner が共有する **幾何・添字・符号規約**を固定する。
**全 runner はこの規約に従って `metric_frame` と `structure_constants` を定義すること。**

⇒ [English version](CONVENTIONS.md) | [SymPy ガイドライン](SymPy_guideline_ja.md)

## 1. 作用域と前提

* ここで扱う量はすべて **フレーム（正規直交基底）** 上の成分で表す。
* `metric_frame` はフレーム計量 $(g_{ab})$ であり、ユークリッド経路では $(g_{ab}=\delta_{ab})$（`Matrix.eye(dim)`）、Lorentzian/LZ-native 経路では $\eta_{ab}=\mathrm{diag}(-1,+1,+1,+1)$。
* `BaseFrameEngine` の通常パイプラインは、**フレーム方向微分が不要になる状況**を主対象にしている。
  具体的には、構造定数 $(C^{a}{}\_{bc})$ と接続係数 $(\Gamma^{a}{}\_{bc})$ が「フレームに関して定数扱い」になる（左不変フレームなど）設定を runner が採用する。
  * **Note:** 低レベル関数 `compute_riemann_tensor(..., frame_deriv=...)` は 局所 spin-2 jet 解析で使うフレーム方向微分を受け取れる。通常パイプラインではこの引数を渡さない。

## 2. 添字・配列のインデックス順

配列格納は以下で固定する：

* 構造定数：`C[a,b,c] = C^a_{bc}`
* 接続係数：`Gamma[a,b,c] = Γ^a_{bc}`
* リーマン曲率：`Riemann[a,b,c,d] = R^a_{bcd}`

添字の意味：

* $(a)$：出力（上付き）成分
* $(b,c,d)$：入力（下付き）成分
  とくに $(\Gamma^{a}{}\_{bc})$ は $(\nabla_{E_c} E_b = \Gamma^{a}{}\_{bc} E_a)$ に対応する（最後の $(c)$ が “微分方向”）。

## 2.5. Lorentzian/LZ-native 署名規約

Lorentzian 縮約・局所 jet 解析では、次の設定を標準とする：

| 項目 | 規約 |
|---|---|
| `signature` | `"lorentzian"` |
| `frame_convention` | `"lz_native"` |
| フレーム計量 | $\eta_{ab}=\mathrm{diag}(-1,+1,+1,+1)$ |
| 時間添字 | `0` |
| 空間添字 | `{1, 2, 3}` |

この経路では、添字の上げ下げを恒等操作として扱ってはならない。`contortion`、Weyl scalar、Pontryagin/Hodge など、符号に敏感な量は必ず `metric` / `metric_inv` を通じて contraction する。

実装上の対応：

* `make_config(..., signature="lorentzian")` は既定で `frame_convention="lz_native"` を選ぶ。
* `dppu/topology/lz_adapter.py` は legacy 配列から LZ-native 添字順への変換を担当する。
* Lorentzian Pontryagin は `dppu/curvature/pontryagin_lz.py` の metric-aware 実装を使い、`epsilon_4d` ではなく `epsilon_tensor_up/down` を使う。
* `dppu/curvature/hodge.py` の `hodge_dual_2form` は 2-form の Hodge 双対を metric-aware に計算する。

## 3. フレーム・コフレームと構造定数の定義（ここが最重要）

フレームの双対を $(\{E_a\})$、コフレーム（1-forms）を $(\{e^a\})$ とする。

### 3.1 コフレームの構造方程式（固定）

$$
de^a = \frac12 C^a{}_{bc} e^b\wedge e^c,
\qquad C^a{}_{bc} = - C^a{}_{cb}.
$$

### 3.2 双対フレームの交換関係（同値）

上の定義は次と同値：

$$
[E_b, E_c] = - C^{a}{}\_{bc} E_a.
$$

> 注意：多くの教科書では $([E_b,E_c]=+f^{a}{}\_{bc}E_a)$ を採用する。
> 本プロジェクトでは **その $(f^{a}{}\_{bc})$ に対して $(C^{a}{}\_{bc}=-f^{a}{}\_{bc})$** の規約を採用している。

### 3.3 runner 実装ルール（推奨）

* **C は手打ちしない**。可能なら runner 側で $de^a$ を明示し、係数比較で $C^a_{bc}$ を抽出して `self.data['structure_constants']=C` に入れる。
* 最低限、 $C^{a}_{bc}$ が **b,c で反対称**になっていることを自動チェックすること。

## 4. 接続（スピン接続）とメトリック整合

接続 1-form を

$$
\omega^{a}{}\_{b} = \Gamma^{a}{}\_{bc} e^c
$$

で定義する。

メトリック整合（ローレンツ接続／直交接続）を仕様として固定：

$$
\omega_{ab} = -\omega_{ba}
\quad(\Leftrightarrow\quad
\Gamma_{abc} = -\Gamma_{bac})
$$

ただし $(\Gamma_{abc} = g_{ad}\Gamma^{d}{}\_{bc})$。

## 5. Levi-Civita 接続（現在の engine の一般 Koszul 実装）

フレームが正規直交で、上記の構造定数規約を採用したとき、Levi-Civita 接続は現在の engine では次の**一般 Koszul 公式**で計算する：

$$
\Gamma^a{}_{bc}
= \frac12\Big(
C^a{}_{bc} + C^c{}_{ba} - C^b{}_{ac}
\Big).
$$

（これは本書 3.2 の交換関係の符号を採用した場合の形である。）

**重要な注記:**

1. この公式は **bi-invariant 計量を仮定しない**。
   左不変フレーム上の Levi-Civita 接続として、Nil³ のような非 bi-invariant な場合でも正しく機能する。

2. SU(2) のように（低い添字で）構造定数が全反対称 $C_{abc} = -C_{bac} = -C_{acb}$ になる特殊な場合は、
   この公式は $\Gamma^a_{bc} = \frac{1}{2} C^a_{bc}$ に簡約される。

3. engine は計算後に **metric compatibility** $\Gamma_{abc} + \Gamma_{bac} = 0$ を自動検証する。
   これに違反する場合は実装エラーとして即座に例外を投げる。

## 6. ねじれと曲率

ねじれ 2-form：

$$
T^a = de^a + \omega^{a}{}\_b \wedge e^b,
\qquad
T^a = \frac12 T^{a}{}\_{bc} e^b\wedge e^c.
$$

曲率 2-form：

$$
R^a{}_b = d\omega^{a}{}\_b + \omega^{a}{}\_c\wedge \omega^{c}{}\_b,
\qquad
R^a{}_b = \frac12 R^{a}{}\_{bcd} e^c\wedge e^d.
$$

## 7. 曲率成分の計算式（engine が実際に使う形）

通常パイプラインの $R^{a}_{bcd}$ は次の形を用いる：

$$
R^{a}{}\_{bcd} = \Gamma^{a}{}\_{ec}\Gamma^{e}{}\_{bd} -\Gamma^{a}{}\_{ed}\Gamma^{e}{}\_{bc} +\Gamma^{a}{}\_{be} C^{e}{}\_{cd}.
$$

> 重要：一般にはここにフレーム方向微分項
> $(E_c(\Gamma^{a}{}\_{bd}) - E_d(\Gamma^{a}{}\_{bc}))$
> が現れる。`BaseFrameEngine` の通常パイプラインではこの項を渡さないため、runner は左不変フレーム等により **$(\Gamma)$ がフレーム方向で定数扱い**になる設定を採用する。
> 局所 jet スクリプトのように時間方向微分が必要な場合は、低レベル関数 `compute_riemann_tensor(..., frame_deriv=...)` を直接用いる。

## 8. 必須セルフチェック（runner が満たすべき整合性）

runner は以下を満たすこと（落ちたら定義が engine と不整合）：

1. 構造定数の反対称：

$$
C^{a}{}\_{bc} + C^{a}{}\_{cb} = 0.
$$

2. メトリック整合（直交接続）：

$$
\omega_{ab} + \omega_{ba} = 0.
$$

3. リーマンの反対称（engine の strict check 対象）：

$$
R_{ab cd} = -R_{ba cd},\qquad
R_{ab cd} = -R_{ab dc}.
$$

---

## 9. Weyl テンソルと共形スカラー

### 9.1 Weyl テンソルの定義（4次元）

$$
C_{abcd} = R_{abcd} - \frac{1}{2}(g_{ac}R_{bd} - g_{ad}R_{bc} - g_{bc}R_{ad} + g_{bd}R_{ac}) + \frac{R}{6}(g_{ac}g_{bd} - g_{ad}g_{bc}).
$$

主要な性質：
- **無跡**： $C^{a}{}\_{bad} = 0$（全添字対で成立）
- **共形不変**：計量の共形変換 $g_{ab} \to \Omega^2 g_{ab}$ のもとで不変
- **共形平坦判定**： $C_{abcd} = 0 \Leftrightarrow$ 共形的平坦

フレーム基底（正規直交）での注意： $g_{ab} = \delta_{ab}$ のため添字の上げ下げは恒等変換に等しく、
 $C_{abcd} = C^{abcd}$ として直接成分を扱える。

### 9.2 Weyl スカラー

$$
C^2 = C_{abcd}\,C^{abcd} = \sum_{a,b,c,d} C_{abcd}^2.
$$

フレーム基底での高速計算：添字上げ下げを省略し、直接2乗和をとる。
$C^2 = 0 \Leftrightarrow$ 共形的平坦（等方的 $S^3 \times S^1$ で $\varepsilon = 0$ のとき成立）。

### 9.3 engine での実装位置

- モジュール：`dppu/curvature/weyl.py`
  - `compute_weyl_tensor(R_abcd, Ricci, R_scalar, metric, dim)` → $C_{abcd}$
  - `compute_weyl_scalar(C_abcd, metric_inv, dim)` → $C^2$
- パイプラインステップ：`E4.3b`（Levi-Civita 曲率計算の直後）

---

## 10. Squashed 等質空間と $\varepsilon$-パラメータ

### 10.1 Squashing の定義

体積保存の異方性変形パラメータ $\varepsilon$ を導入する。
$\varepsilon = 0$ が等方（アイソトロピック）基準点。物理的範囲： $\varepsilon \in (-1, +\infty)$（ $\varepsilon = -1$ で構造定数が発散し特異点）。

### 10.2 トポロジー別の構造定数スケーリング

左不変フレームの基底構造定数 $C^a{}_{bc}(\varepsilon=0)$ を以下のスケーリングで変形する。

**$S^3 \times S^1$（SU(2)）:**

| フレーム添字 $a$ | スケール因子 $\lambda_a(\varepsilon)$ |
|---|---|
| $a \in \{0, 1\}$ | $(1+\varepsilon)^{2/3}$ |
| $a = 2$ | $(1+\varepsilon)^{-4/3}$ |

体積保存の確認： $\lambda_0\lambda_1\lambda_2 = (1+\varepsilon)^{2/3+2/3-4/3} = 1$。

**$\mathrm{Nil}^3 \times S^1$（Heisenberg 群）:**

| 非自明な添字 $a$ | スケール因子 |
|---|---|
| $a = 2$（非可換成分） | $(1+\varepsilon)^{-4/3}$ |

$\varepsilon \to +\infty$ で構造定数が消失し、平坦 $T^3$ 的挙動に漸近する。

**$T^3 \times S^1$（Abelian 群）:**

構造定数は恒等的にゼロ（ $C^a{}_{bc} = 0$）。 $\varepsilon$ による変形は定義されない。
$C^2 = 0$ が全域で成立し、Weyl 項は常にゼロとなる（null test）。

### 10.3 物理的極限

| 極限 | 物理的意味 |
|---|---|
| $\varepsilon = 0$ | 等方的 $S^3$ ： $C^2 = 0$、Paper I の安定真空 |
| $\varepsilon \to +\infty$ | $\mathrm{Nil}^3$ の平坦極限（ $C^2 \to 0$ 漸近） |
| $\varepsilon \to -1^+$ | $S^3$ 構造の特異点（物理的到達不能） |
| $\varepsilon = -2$ | $C^2 = 0$ の数学的根だが $\varepsilon < -1$ のため物理的除外 |

> $\varepsilon$ は完全な spin-2 五重項における $T_1$ 方向の座標に対応する（§11 参照）。相転移は $\varepsilon_{c+} \approx 0.483$（squash softening）および $\varepsilon_{c-} \approx -0.293$（shear 凝縮）で起こる。

---

## 11. 完全 shear テンソル基底と Spin-2 五重項

### 11.1 traceless 対称テンソル基底 $T_A$

$S^3$ の完全 shear 変形は traceless 対称 $3\times3$ 行列 $h_{ij}$（spin-2 五重項）で表される。正規直交基底 $\{T_A\}\_{A=1}^{5}$（ $\mathrm{Tr}(T_A T_B) = \delta\_{AB}$）：

| 添字 $A$ | 行列 | SO(2) チャージ（squash 軸 = 2 軸） |
|---|---|---|
| 1 | $\mathrm{diag}(1,1,-2)/\sqrt{6}$ | $m=0$（singlet） |
| 2 | $\mathrm{diag}(1,-1,0)/\sqrt{2}$ | $m=\pm 2$ doublet |
| 3 | $(e_0 e_1^T + e_1 e_0^T)/\sqrt{2}$ | $m=\pm 2$ doublet |
| 4 | $(e_0 e_2^T + e_2 e_0^T)/\sqrt{2}$ | $m=\pm 1$ doublet |
| 5 | $(e_1 e_2^T + e_2 e_1^T)/\sqrt{2}$ | $m=\pm 1$ doublet |

コード例（SymPy）：
```python
T[1] = Matrix([[1, 0, 0], [0, 1, 0], [0, 0, -2]]) / sqrt(6)
T[2] = Matrix([[1, 0, 0], [0, -1, 0], [0, 0, 0]]) / sqrt(2)
T[3] = Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 0]]) / sqrt(2)
T[4] = Matrix([[0, 0, 1], [0, 0, 0], [1, 0, 0]]) / sqrt(2)
T[5] = Matrix([[0, 0, 0], [0, 0, 1], [0, 1, 0]]) / sqrt(2)
```

### 11.2 $(\varepsilon, s)$ パラメータの五重項への埋め込み

squash パラメータ $\varepsilon$ と shear パラメータ $s$ は、それぞれ $T_1$、 $T_2$ 方向の座標に対応する：

$$
h = q_1 T_1 + q_2 T_2 + \cdots,\quad
q_1 = \tfrac{\sqrt{6}}{2}\ln(1+\varepsilon),\quad
q_2 = \sqrt{2}\ln(1+s).
$$

`S3Geometry` が用いる 2 パラメータ・フレームの scale factor：
```
factor_0 = (1+ε)^{2/3} (1+s)²
factor_1 = (1+ε)^{2/3} / (1+s)²
factor_2 = (1+ε)^{−4/3}          （s に依存しない）
```

### 11.3 等方点と Schur 縮退

等方点（ $\varepsilon = s = 0$）では SO(3) 対称性により Hessian が単位行列に比例する（Schur の補題）：

$$
H_{AB} = \mu_q^2\,\delta_{AB},\quad \mu_q^2 = 48\pi^2.
$$

等方点での off-diagonal Hessian $H_{AB}$（ $A \neq B$）および spin-2 と他 sector（spin-0, spin-1）間の cross-Hessian は全てゼロ。

### 11.4 SO(3) → SO(2) 分裂と相転移

$\varepsilon \neq 0$ では対称性が SO(2) に下がり、5 重縮退が最大 3 固有値に分裂する：

| Sector | 臨界点 | 物理的意味 |
|---|---|---|
| $m=0$ singlet（ $T_1$） | $\varepsilon_{c+} \approx 0.483$ | squash softening（第 1 soft mode） |
| $m=\pm2$ doublet（ $T_2, T_3$） | $\varepsilon_{c-} \approx -0.293$ | $s \neq 0$ への 2 次相転移 |
| $m=\pm1$ doublet（ $T_4, T_5$） | 両 $\varepsilon_{c\pm}$ で massive | 相転移なし |

---

## 12. Spin-1 セクター：Physical/Auxiliary 分解

### 12.1 生変数

Spin-1 セクターは 6 個の生 DOF を持つ：twist $\omega_k$ と mixing $\delta_k$（ $k=0,1,2$）。エンジン内では `omega1, omega2, omega3` および `delta0, delta1, delta2`（BOTH fiber mode）として定義される。

### 12.2 null mode と non-dynamical auxiliary

field-space metric $G$（運動エネルギー計量、rank-3 / 6D）の null 方向：

$$
Y_k \propto \left(\omega_k,\, \delta_k\right) = \left(1,\, \frac{2L}{r_0}\right).
$$

**分類 — Case B（確認済み）**：$Y_k$ は non-dynamical auxiliary（一次拘束と類似の）方向であり、ゲージ対称性でも変数の冗長性でもない（ $Y_k$ 方向に沿って $\delta g_{\mu\nu} \neq 0$）。

### 12.3 物理的三重項 $X_k$

伝播する唯一の spin-1 三重項：

$$
X_k = \frac{-2\omega_k + 3\delta_k}{\sqrt{13}},
$$

一般化固有値 $\lambda_{\rm phys} \approx 4.19 > 0$（3 重縮退、SO(3) 対称）。

### 12.4 一般化固有値問題

Spin-1 セクターの安定性は次式で評価する：

$$
H\,v = \lambda\,G\,v.
$$

手順：
1. $H$（ $V_{\rm eff}$ の $(\omega_k, \delta_k)$ に関する 6×6 Hessian）を有限差分で計算。
2. $G$（6×6 field-space metric）を velocity-mode engine（`enable_velocity=True`）で計算。
3. $G$ の null 方向 $Y_k$ を同定し、直交補空間に射影。
4. 3D 物理部分空間で $H_{\rm phys} v = \lambda_{\rm phys} G_{\rm phys} v$ を解く。

> **注意**： $H$ の 6×6 問題で出る負固有値は $G$ の null 方向に乗っており、伝播しない。タキオンではない。

---

## 13. UnifiedEngine と DOFConfig アーキテクチャ

### 13.1 クラス階層

現在のエンジンは全トポロジーを統一インターフェースに集約している：

```
BaseFrameEngine          (dppu/engine/pipeline.py)
    └── TopologyEngine   (dppu/topology/base_topology.py)  ← abstract
            ├── S3Geometry    (dppu/topology/s3.py)
            ├── T3Geometry    (dppu/topology/t3.py)
            ├── Nil3Geometry  (dppu/topology/nil3.py)
            └── Sol3Geometry  (dppu/topology/sol3.py)
```

`UnifiedEngine`（`dppu/topology/unified.py`）は `DOFConfig` に基づいて適切なサブクラスを選択するファクトリである。

### 13.2 DOFConfig パラメータ

| パラメータ | 型 | 意味 |
|---|---|---|
| `topology` | `TopologyType` | `S3` / `T3` / `NIL3` / `SOL3` |
| `enable_squash` | `bool` | $\varepsilon$ 変形を有効化 |
| `enable_shear` | `bool` | $s$（ $T_2$）shear を有効化 |
| `enable_offdiag_shear` | `bool` | $q_3,q_4,q_5$ off-diagonal shear を有効化 |
| `fiber_mode` | `FiberMode` | `NONE` / `TWIST` / `MIXING` / `BOTH` |
| `isotropic_twist` | `bool` | 3 方向で共通の $\omega$ を使用 |
| `torsion_mode` | `Mode` | `AX` / `VT` / `MX` |
| `ny_variant` | `NyVariant` | Nieh-Yan バリアント |
| `enable_velocity` | `bool` | $G$ 計量計算用の velocity シンボルを有効化 |
| `skip_antisymmetry_check` | `bool` | 重い Riemann 反対称性検証を省略 |
| `weyl_source` | `CurvatureSource` | action の Weyl 項に `LC` / `EC` のどちらを使うか |
| `riemann_source` | `CurvatureSource` | lambdify など診断用 Riemann の既定ソース |
| `signature` | `str` | `"euclidean"` / `"lorentzian"` |
| `frame_convention` | `str` | `"legacy_euclidean"` / `"lz_native"` |

### 13.3 使用例

```python
from dppu.topology import UnifiedEngine, TopologyType, FiberMode, make_config
from dppu.torsion.mode import Mode
from dppu.torsion.nieh_yan import NyVariant

cfg = make_config(
    TopologyType.S3,
    enable_squash=False,
    fiber_mode=FiberMode.BOTH,
    isotropic_twist=False,
    torsion_mode=Mode.MX,
    ny_variant=NyVariant.FULL,
    signature="lorentzian",
)
engine = UnifiedEngine(cfg)
engine.run()

Veff = engine.data['potential']
fp   = engine.get_free_params()   # 有効な SymPy Symbol の dict
```

### 13.4 レガシーエンジンとの対応

| レガシークラス | DOFConfig 等価設定 |
|---|---|
| `S3S1Engine` | `S3, squash=True, shear=False, fiber=NONE` |
| `T3S1Engine` | `T3, squash=False, shear=False, fiber=NONE` |
| `Nil3S1Engine` | `NIL3, squash=True, shear=False, fiber=NONE` |

現在の実装ではレガシークラス名からの自動生成 API は提供しない。新規コードは `make_config(...)` または `make_engine(...)` で明示的に設定する。

