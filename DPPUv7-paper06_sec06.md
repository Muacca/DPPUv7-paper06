## 6. Four-topology admissibility classification

Section 3 の Hamiltonian admissibility と Section 4 の P-channel diagnostics を統合した local admissibility classification を Table 6 にまとめる。

**Table 6.** Four-topology admissibility classification.

| topology | EH | AX | VT | MX |
|----------|----|----|----|----|
| $S^3$ | `L_ADMISSIBLE` | `L_ADMISSIBLE` | `L_ADMISSIBLE` | `L_CONDITIONALLY_ADMISSIBLE` |
| $T^3$ | `L_ADMISSIBLE` | `L_ADMISSIBLE` | `L_ADMISSIBLE` | `L_CONDITIONALLY_ADMISSIBLE` |
| $Nil^3$ | `L_ADMISSIBLE` | `L_ADMISSIBLE` | `L_ADMISSIBLE` | `L_CONDITIONALLY_ADMISSIBLE` |
| $Sol^3$ | `L_ADMISSIBLE` | `L_ADMISSIBLE` | `L_ADMISSIBLE` | `L_CONDITIONALLY_ADMISSIBLE` |

`L_ADMISSIBLE` の label は、Hamiltonian constraint の closure（Section 3.2）、active torsion の auxiliary closure（`AUX_PASS`）、および $P_{\rm form}=0$, $P_{\rm int}=0$ の diagnostic check を通過したことを表す。EH は pure curvature reference として P-channel が適用されない。`L_CONDITIONALLY_ADMISSIBLE` の label は、MX branch が auxiliary real branch 条件（ $\kappa^4\theta_{\rm NY}^2+1>0$）と auxiliary shell 解釈（ $\eta=V=0$）を前提として初めて admissible となることを意味する。off-shell の $P_{\rm int}^{\rm MX}$ および $Q_{\rm NY}$ は diagnostic / boundary channel として記録されるが、admissibility の obstruction ではない。

この表の構造的要点は三点に要約できる。第一に、EH/AX/VT/MX の mode pattern は四つの topology にわたって一致する。第二に、MX の conditional admissibility の条件内容も四つの topology に共通であり（Section 3.2 の auxiliary 行列式形参照）、topology は条件の有無を変えない。第三に、 $Sol^3$ は compact quotient に関する caveat を持つが（Appendix A 参照）、main classification table では他の三つの topology と同じ行構造を持つ。

この topology-robust pattern は Section 5 の χ-universality の背景をなす。admissibility 自体が topology によって変化しないため、topology dependence は admissibility label ではなく、diagnostic/orbit 側の量（ $C_{\rm topology}$ および $\chi$）に現れる。この意味で、Section 5 の bridge result は admissibility classification の上に自然に乗る帰結である。
