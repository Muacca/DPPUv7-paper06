## 3. Hamiltonian admissibility by topology

### 3.1 Admissibility criterion

本稿では、reduced sector における admissibility を次の基準で判定する。

1. **EH constraint closure**: pure EH 部分の Hamiltonian constraint が Lorentzian convention のもとで閉じること（`EH_PASS`）。
2. **Reduced action skeleton**: 各 torsion mode の EC+NY reduced Lagrangian が well-defined な velocity structure を持つこと（`LRED_PASS`）。
3. **Torsion auxiliary closure**: auxiliary field $\eta$, $V$ の algebraic equations が実数解を持つこと（`AUX_PASS` または `AUX_CONDITIONAL_REAL_BRANCH`）。

この admissibility は、orbit stability、unitarity、または full constraint closure の主張ではない。

### 3.2 Hamiltonian constraint

lapse variation $\delta S / \delta N = 0$ より Hamiltonian constraint

$$
\mathcal{H} \equiv A_\Sigma\, q\!\left(F_\Sigma + \frac{\dot q^2}{N^2}\right) = 0,\qquad
\text{すなわち}\quad F_\Sigma + \frac{\dot q^2}{N^2} = 0
$$

を得る。EH branch（ $\eta=V=0$）では $F_\Sigma\big|_{\rm EH} = \chi$（各 topology で Table 1 の $\chi$ を代入）となり、これが `EH_PASS` の要件である。具体的には：

| topology | $F_{\rm EH}$ (on-shell form) | constraint residual |
|----------|------------------------------|---------------------|
| $S^3$ | $4$ | $4N^2 + \dot r^2 = 0$ |
| $T^3$ | $0$ | $\dot r^2 = 0$ |
| $Nil^3$ | $-1/12$ | $N^2 - 12\dot R^2 = 0$ |
| $Sol^3$ | $-1/3$ | $3N^2 - 9\dot q^2 = 0$ (local/iso) |

torsion auxiliary equations は $\delta S/\delta\eta=0$, $\delta S/\delta V=0$ から得られる代数方程式である。AX では $\eta$ のみ、VT では $V$ のみが非自明に現れ、それぞれ `AUX_PASS` となる。MX では $\eta$, $V$ の連立 auxiliary 行列式が

$$
\det\!\begin{pmatrix}\partial_\eta(\partial_\eta F) & \partial_V(\partial_\eta F)\\ \partial_\eta(\partial_V F) & \partial_V(\partial_V F)\end{pmatrix}
= -\frac{4q^2(\kappa^4\theta_{\rm NY}^2+1)}{9} < 0
$$

となる。 $q>0$ の reduced-sector domain と real $\kappa,\theta_{\rm NY}$ のもとでは $\kappa^4\theta_{\rm NY}^2+1>0$ なので auxiliary Hessian は非退化である。本稿では、この mixed auxiliary channel を real branch の選択と origin tagging を伴う branch として保持するため、`AUX_CONDITIONAL_REAL_BRANCH`、すなわち MX の `H_CONDITIONAL` として記録する。
この負の determinant は安定な minimum を意味しない。ここで用いる `conditional admissibility` は、non-degenerate mixed auxiliary saddle における algebraic closure と real branch selection のみを要求する弱い基準であり、Hessian の符号（dynamical stability）とは独立である。

### 3.3 Mode status

以上の判定をまとめると Table 2 を得る。結果は四つの topology に対して一致する。

**Table 2.** Hamiltonian status by topology and mode.

| mode | Hamiltonian status | interpretation |
|------|--------------------|----------------|
| EH | `H_PASS` | pure EH reference branch |
| AX | `H_PASS` | active torsion admissible |
| VT | `H_PASS` | active torsion admissible |
| MX | `H_CONDITIONAL` | real auxiliary branch を経由した conditional admissibility |

MX の `H_CONDITIONAL` は obstruction や failure ではない。real auxiliary branch が存在する前提のもとで、後続の P-channel diagnostics および $Q_{\rm NY}$ boundary channel の off-shell 振る舞いを diagnostic として記録することが、admissibility の条件内容である。

**注.** `H_PASS` および `L_ADMISSIBLE` は dynamical stability を意味しない。 $Sol^3$ の table entry は local/isotropic reduction の範囲内のものであり、compact quotient の完全解決を意味しない。
