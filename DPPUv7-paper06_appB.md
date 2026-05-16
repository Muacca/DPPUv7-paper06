## Appendix B. Computational checks and reproducibility

### B.1 Algebraic core of the $\chi$-bridge identity

Section 5.2 の family-level identity $C+9\chi=0$ の代数的核心を記す。補助量 $\Delta$ を

$$
\Delta = a^2-2ab-2ac+b^2-2bc+c^2+4u^2+4uv+4v^2
$$

と定義すると、 $\chi$ および $C$ の二つの量は

$$
\chi=-\frac{\Delta}{12},\qquad C=\frac{3\Delta}{4}
$$

と表される。よって

$$
C+9\chi=\frac{3\Delta}{4}-\frac{9\Delta}{12}=0.
$$

この恒等式は $\Delta$ の任意の値に対して成立し、特定の topology を仮定しない。四つの topology の特殊化値（Table A1）はその代入例である。

### B.2 Full topology branch taxonomy

全 16 branch（ $S^3/T^3/Nil^3/Sol^3$ $\times$ EH/AX/VT/MX）の 主な admissibility classification を Table B1 に示す。Table 6（Section 6）はこの分類 pattern の概要である。
完全な $F$-branch 式、Hamiltonian constraint、auxiliary equations、auxiliary shell 解は `admissibility_classification.py` にて再現可能である。

**Table B1.** Compact full branch taxonomy.

| topology | mode | $\chi$ | $F_{\rm shell}$ | Hessian determinant | label |
|----------|------|--------|-----------------|---------------------|-------|
| $S^3$ | EH | $4$ | $4$ | n/a | L_ADMISSIBLE |
| $S^3$ | AX | $4$ | $4$ | $-2$ | L_ADMISSIBLE |
| $S^3$ | VT | $4$ | $4$ | $2q^2/9$ | L_ADMISSIBLE |
| $S^3$ | MX | $4$ | $4$ | $-4q^2(\kappa^4\theta_{\rm NY}^2+1)/9$ | L_CONDITIONALLY_ADMISSIBLE |
| $T^3$ | EH | $0$ | $0$ | n/a | L_ADMISSIBLE |
| $T^3$ | AX | $0$ | $0$ | $-2$ | L_ADMISSIBLE |
| $T^3$ | VT | $0$ | $0$ | $2q^2/9$ | L_ADMISSIBLE |
| $T^3$ | MX | $0$ | $0$ | $-4q^2(\kappa^4\theta_{\rm NY}^2+1)/9$ | L_CONDITIONALLY_ADMISSIBLE |
| $Nil^3$ | EH | $-1/12$ | $-1/12$ | n/a | L_ADMISSIBLE |
| $Nil^3$ | AX | $-1/12$ | $-1/12$ | $-2$ | L_ADMISSIBLE |
| $Nil^3$ | VT | $-1/12$ | $-1/12$ | $2q^2/9$ | L_ADMISSIBLE |
| $Nil^3$ | MX | $-1/12$ | $-1/12$ | $-4q^2(\kappa^4\theta_{\rm NY}^2+1)/9$ | L_CONDITIONALLY_ADMISSIBLE |
| $Sol^3$ | EH | $-1/3$ | $-1/3$ | n/a | L_ADMISSIBLE |
| $Sol^3$ | AX | $-1/3$ | $-1/3$ | $-2$ | L_ADMISSIBLE |
| $Sol^3$ | VT | $-1/3$ | $-1/3$ | $2q^2/9$ | L_ADMISSIBLE |
| $Sol^3$ | MX | $-1/3$ | $-1/3$ | $-4q^2(\kappa^4\theta_{\rm NY}^2+1)/9$ | L_CONDITIONALLY_ADMISSIBLE |

### B.3 Full orbit-sheet atlas

Reduced vacuum orbit atlas の 24 branch/sheet rows を Table B2 に示す。 $S^3$ と $T^3$ は各 branch 一行、 $Nil^3$ と $Sol^3$ は expanding/collapsing の二 sheets を持つ。Table 8（Section 7.3）はこの branch/sheet count を要約したものである。
各 row の orbit equation、constraint residual、数値 monitor は `orbit_atlas.py` にて再現可能である。

**Table B2.** Compact full orbit-sheet atlas.

| topology | mode | sheet | $\dot q_0$ | orbit class | check |
|----------|------|-------|------------|-------------|-------|
| $S^3$ | EH | none | no real solution | NO_REAL_AUX_SHELL_ORBIT | PASS |
| $S^3$ | AX | none | no real solution | NO_REAL_AUX_SHELL_ORBIT | PASS |
| $S^3$ | VT | none | no real solution | NO_REAL_AUX_SHELL_ORBIT | PASS |
| $S^3$ | MX | none | no real solution | NO_REAL_AUX_SHELL_ORBIT | PASS |
| $T^3$ | EH | zero | $0$ | STATIC_OR_DEGENERATE | PASS |
| $T^3$ | AX | zero | $0$ | STATIC_OR_DEGENERATE | PASS |
| $T^3$ | VT | zero | $0$ | STATIC_OR_DEGENERATE | PASS |
| $T^3$ | MX | zero | $0$ | STATIC_OR_DEGENERATE | PASS |
| $Nil^3$ | EH | expanding | $\sqrt{3}/6$ | MONOTONIC_EXPANSION | PASS |
| $Nil^3$ | EH | collapsing | $-\sqrt{3}/6$ | SINGULAR_APPROACH | PASS |
| $Nil^3$ | AX | expanding | $\sqrt{3}/6$ | MONOTONIC_EXPANSION | PASS |
| $Nil^3$ | AX | collapsing | $-\sqrt{3}/6$ | SINGULAR_APPROACH | PASS |
| $Nil^3$ | VT | expanding | $\sqrt{3}/6$ | MONOTONIC_EXPANSION | PASS |
| $Nil^3$ | VT | collapsing | $-\sqrt{3}/6$ | SINGULAR_APPROACH | PASS |
| $Nil^3$ | MX | expanding | $\sqrt{3}/6$ | MONOTONIC_EXPANSION | PASS |
| $Nil^3$ | MX | collapsing | $-\sqrt{3}/6$ | SINGULAR_APPROACH | PASS |
| $Sol^3$ | EH | expanding | $\sqrt{3}/3$ | MONOTONIC_EXPANSION | PASS |
| $Sol^3$ | EH | collapsing | $-\sqrt{3}/3$ | SINGULAR_APPROACH | PASS |
| $Sol^3$ | AX | expanding | $\sqrt{3}/3$ | MONOTONIC_EXPANSION | PASS |
| $Sol^3$ | AX | collapsing | $-\sqrt{3}/3$ | SINGULAR_APPROACH | PASS |
| $Sol^3$ | VT | expanding | $\sqrt{3}/3$ | MONOTONIC_EXPANSION | PASS |
| $Sol^3$ | VT | collapsing | $-\sqrt{3}/3$ | SINGULAR_APPROACH | PASS |
| $Sol^3$ | MX | expanding | $\sqrt{3}/3$ | MONOTONIC_EXPANSION | PASS |
| $Sol^3$ | MX | collapsing | $-\sqrt{3}/3$ | SINGULAR_APPROACH | PASS |

### B.4 Summary of main computational checks

Table B3 に logical verification target を示す。最初の二つの check は同一の symbolic script で評価されるが、これは family-level identity とその topology specialization が同じ導出済み $(a,b,c,u,v)$ data を共有しているためである。

**Table B3.** 計算 check とその役割。

| check | content | related section | script |
|-------|---------|------|--------|
| $\chi$-bridge derivation | $S^3/T^3/Nil^3/Sol^3$ の $(a,b,c,u,v)$ を抽出 、LC spatial curvature から $\chi$ を導出し、EC+NY MX curvature から raw $P_{\rm int}^{\rm MX}$ を計算して $C$ を抽出 | Section  4.3, 5.2 | `chi_bridge_symbolic.py` |
| $P_{\rm form}$ cancellation | AX/VT/MX の各 torsion branch で EC+NY engine から curvature を構成し $P_{\rm form}$ を計算 | Section 4.2 | `pform_cancellation.py` |
| admissibility classification | engine-derived $F_{\rm branch}$ から lapse constraint と auxiliary closure を symbolic に導出 | Section 3, 6 | `admissibility_classification.py` |
| orbit atlas | engine-derived $F_{\rm branch}$ の auxiliary shell から orbit equation と orbit class を導出 | Section 7 | `orbit_atlas.py` |
