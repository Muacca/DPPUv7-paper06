## 4. P-channel dictionary and topology response

### 4.1 Channel separation

本稿では Pontryagin 型の量を次の三つに分離する。

**Table 3.** P-channel dictionary.

| channel | result | role | distinction |
|---------|--------|------|-------------|
| $P_{\rm form}$ | active torsion branches で exact zero | structural cancellation | dynamical / orbit stability とは別 |
| $P_{\rm int}$ | MX internal-pair diagnostic；$C_{\rm topology}=-9\chi$ | diagnostic topology response | true Pontryagin density ではない |
| $Q_{\rm NY}$ | endpoint / boundary channel；VT zero, AX on-shell trivial, MX boundary-conditional | boundary bookkeeping | bulk の動力学的機構ではない |

$P_{\rm int}$ は orthonormal frame の internal index pair の contraction に基づく diagnostic quantity であり、form-Hodge Pontryagin density [13] ではない。この区別は本稿の解釈において本質的である。

### 4.2 $P_{\rm form}$ cancellation by block orthogonality

$P_{\rm form}$ は curvature 2-form の form-Hodge contraction として定義する。成分表示では

$$
P_{\rm form}=
\frac{1}{4}\sum_{a,b,c,d}R^{ab\,cd}\,(*R_{ab})_{cd},
\qquad
(*R_{ab})_{cd}=
\frac{1}{2}\epsilon_{cd\,ef}R_{ab}{}^{ef}
$$

である。ここで Hodge dual は Lorentzian metric-aware な form-index dual として取る。

**Lemma 1 (block-orthogonality cancellation).** 本稿の Lorentzian isotropic EC+NY reduced sector では、AX/VT/MX の各 active torsion branch および $S^3/T^3/Nil^3/Sol^3$ の各 topology に対して、form-Hodge Pontryagin density は恒等的に

$$
P_{\rm form}=0
$$

となる。

**Proof sketch.** reduced coframe $e^0=N\,dt$, $e^i=q\,\sigma^i$ のもとで、非零 curvature 2-form は同一 internal pair ごとに lapse-space block と purely spatial block に分離される。form-Hodge dual は各 block を補完 block に写すが、同一 internal pair で縮約可能な相手が存在しないため、pairwise contraction が消える。全 branch の symbolic verification は Appendix B に記録する。この cancellation は本稿の isotropic reduced setting における構造的な結果であり、full theory や anisotropic extension には一般には成立しない。また、 $P_{\rm form}=0$ は dynamical stability の主張ではない。

### 4.3 $P_{\rm int}^{\rm MX}$ as internal-pair diagnostic

AX および VT では、 $P_{\rm int}$ の off-shell 表式が恒等的に zero となる。MX では off-shell nonzero channel が現れるが、auxiliary shell（ $\eta=V=0$；real branch 条件のもと）上では zero に戻る。この off-shell 振る舞いは Hamiltonian admissibility を破らない。

MX branch の off-shell $P_{\rm int}^{\rm MX}$ は四つの topology にわたって共通の functional template

$$
P_{\rm int}^{\rm MX}=
\frac{2V\eta\!\left(-V^2q^2+9\eta^2+C_{\rm topology}\right)}{9q^3}
$$

に一致する。topology への依存は係数 $C_{\rm topology}$ にのみ局在し、四つの topology でそれぞれ

**Table 4.** $C_{\rm topology} = -9\chi$ bridge table.

| topology | $\chi$ | $C_{\rm topology}$ | $C_{\rm topology}+9\chi$ |
|----------|-------:|-------------------:|-------------------------:|
| $S^3$ | $4$ | $-36$ | $0$ |
| $T^3$ | $0$ | $0$ | $0$ |
| $Nil^3$ | $-1/12$ | $3/4$ | $0$ |
| $Sol^3$ | $-1/3$ | $3$ | $0$ |

が成立する。 $C_{\rm topology}$ は raw $P_{\rm int}^{\rm MX}$ の symbolic 表式から template fitting により抽出したものであり、 $\chi$ を前提とせずに得られる。 $\chi$ は EH/spatial-curvature sector から独立に導出される（Section 2.2 参照）。この二つの量の一致が Section 5 で述べる χ-bridge result の実質的内容である。

### 4.4 $Q_{\rm NY}$ as endpoint channel

Nieh-Yan primitive [8,13] は $Q_{\rm NY}=e^a\wedge T_a$ であり、 $P_{\rm form}$, $P_{\rm int}$ とは独立の exact-form endpoint quantity として扱う。AX/MX branch は axial torsion $\eta$ を活性に持ち、その reduced pullback は

$$
Q_{\rm NY}=\frac{6\eta}{q}\,e^{123},\qquad
\int_\Sigma Q_{\rm NY}=6\mathcal V_\Sigma q^2\eta,
$$

orientation は $Q(t_f)-Q(t_i)$ で固定する。ここで $\mathcal V_\Sigma$ は unit coframe volume の normalization factor であり、torsion variable $V(t)$ とは区別する。VT branch では axial component $\eta=0$ が EOM から従うため、off-shell から $Q_{\rm NY}=0$ である。auxiliary shell 上では、AX は on-shell trivial、MX は boundary-conditional となる。 $Q_{\rm NY}$ は bulk の動力学的機構ではなく、endpoint / boundary bookkeeping として扱うべき量である。
