## Appendix A. Topology conventions, local reductions, and caveats

### A.1 Structure constants and $\chi$ values

本稿で用いる四つの topology の spatial structure constants を five-parameter LZ-native isotropic family の特殊化として与える。family のパラメータは $(a,b,c,u,v)$ であり、共変指標を反対称化した構造定数を

$$
C^1{}_{23}=\frac{a}{q},\quad C^2{}_{31}=\frac{b}{q},\quad C^3{}_{12}=\frac{c}{q},\quad
C^2{}_{12}=\frac{u}{q},\quad C^3{}_{13}=\frac{v}{q}
$$

と定義する。このとき三次元 Ricci scalar は

$$
{}^{(3)}R\,q^2 = -\tfrac{1}{2}\,\Delta,\qquad
\Delta = a^2-2ab-2ac+b^2-2bc+c^2+4u^2+4uv+4v^2
$$

であり、本稿の $\chi$ は

$$
\chi = -\frac{\Delta}{12}
$$

によって定義される。各 topology の parameter values と $\chi$ を Table A1 に示す。

**Table A1.** Topology conventions, $\chi$ values, and reduction status.

| topology | $(a,b,c,u,v)$ | $\Delta$ | $\chi$ | reduction status | caveat |
|----------|--------------|--------:|-------:|-----------------|--------|
| $S^3$ | $(4,4,4,0,0)$ | $-48$ | $4$ | closed isotropic reference | — |
| $T^3$ | $(0,0,0,0,0)$ | $0$ | $0$ | flat reference | — |
| $Nil^3$ | $(0,0,-1,0,0)$ | $1$ | $-1/12$ | isotropic-scale local branch | anisotropic EOM は将来課題 |
| $Sol^3$ | $(0,0,0,1,-1)$ | $4$ | $-1/3$ | isotropic-scale local branch | compact quotient は本稿範囲外 |

**検証.** $S^3$: $\Delta=(16-32-32+16-32+16)=-48$, $\chi=4$. $Nil^3$: $\Delta=1$, $\chi=-1/12$. $Sol^3$: $\Delta=(4-4+4)=4$, $\chi=-1/3$.

### A.2 Distinguished loci in the five-parameter family

Table A2 は、Table A1 で固定した topology order と $(\Delta,\chi)$ convention から、bridge-facing な追加情報だけを示す。 $C$ は Appendix B.1 の $C=3\Delta/4$ から計算し、各行で $C+9\chi=0$ を満たす。

**Table A2.** Bridge-facing structural loci.

| topology | $C=3\Delta/4$ | structural locus | paper04 echo |
|----------|---------------:|-----------------|--------------|
| $S^3$ | $-36$ | isotropic class-A diagonal | triaxial |
| $T^3$ | $0$ | origin | none |
| $Nil^3$ | $3/4$ | single class-A axis | uniaxial |
| $Sol^3$ | $3$ | balanced class-B diagonal pair | biaxial |

paper04 echo 列は interpretive のみである。paper04 の response dictionary との structural resonance を記録するものであり、機構の等価性や共通の定理を主張するものではない。

### A.3 Reduction scope for $Nil^3$ and $Sol^3$

$Nil^3$ は Heisenberg 群の左不変計量から定まる非 Abelian homogeneous geometry [9,10] であり、その sign convention として $C^3{}_{12}=-1/q$ を採用する（ $c=-1$）。isotropic scale reduction では一つの scale parameter $q$ を全方向に共通とする。独立な scale factors を導入した anisotropic EOM は本稿では導出しない。

$Sol^3$ は solvable Lie 群に基づく geometry [9,10] であり、balanced class-B diagonal pair $u=1$, $v=-1$ として five-parameter family に入る。これは本稿の isotropic-scale reduced entry である。 $Sol^3$ の compact realization には mapping-torus 構成が必要であり、lattice normalization の選択・spin structure の分類・spectrum data の確定は本稿の範囲外である。Table 6 および Table 7 の $Sol^3$ entry は、この local/isotropic reduced branch に限定した結果である。

### A.4 Lorentzian EC+NY action conventions

本稿では Lorentzian signature $(-,+,+,+)$ の orthonormal frame を採用し、Nieh-Yan coupling は real $\theta_{\rm NY}$ によってパラメタライズする。 $P_{\rm form}$ の frame block-diagonal 構造は Lorentzian と Euclidean で共通であるが、Hamiltonian constraint の符号規則は Lorentzian convention に従う（Section 3.2 参照）。 $P_{\rm int}$ の definition と $Q_{\rm NY}$ の orientation は paper05 [6] の conventions を引き継ぐ。
$S^3$ specialization の定数項 $C=-36$ は paper03（"Unified Geometric Landau EFT", [3]）の MX-P expression の定数項と一致する。一方、 $V^2q^2$ 項の符号は、paper03 [3] の Euclidean convention と本稿の Lorentzian metric-aware epsilon / legacy temporal-index map の差に由来する。
