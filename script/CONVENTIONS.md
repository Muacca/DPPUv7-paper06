# Geometry, symbols and normalization

The action and all equations are defined in the manuscript. This guide maps them to the computational interface.

| Symbol | Meaning |
| --- | --- |
| `q` | Positive common scale; physical round radius for S3 |
| `n1,n2,n3` | Constant eigenvalues of the real symmetric class-A structure matrix |
| `H` | Proper-time logarithmic scale derivative, Dq/q |
| `HD,etaDD,bDD` | Proper-time derivatives of H, etaD and bD in exterior differentiation |
| `eta,etaD` | Geometric axial amplitude and its proper-time derivative |
| `b,bD` | Boost connection coefficient H−V/3 and its proper-time derivative |
| `N,qd,qdd,Nd,ed,V,Vd` | Lapse and coordinate-time jets in `reduction.json` |
| `volume,kappa,theta` | Reference spatial volume, positive gravitational parameter and constant NY coupling |

The signature is (-,+,+,+), the positive volume form is e0123, and epsilon0123=+1 while epsilon^0123=−1. Differential forms use increasing index tuples. The tensor is all lowered; a stored coefficient of `e^c wedge e^d` for c<d equals R_abcd. The curvature is `d omega + omega wedge omega`, and torsion is geometric `de + omega wedge e`.

The spatial structure equation is $`d\sigma^i=\tfrac12 n_{ij}\epsilon_{jkl}\sigma^k\wedge\sigma^l`$. Matrix n is constant. Its eigenvalues are (0,0,0), (0,0,−1), (0,1,−1) and (2,2,2) for T3, Nil, Sol and S3. A common scale does not imply isotropic spatial Ricci curvature. General symbolic formulas are local class-A identities; the listed spatial integrals use the four fixed compact quotients.

`NY` and `Pontryagin` are the positive volume-form coefficients of complete dQ and the all-internal-index sum R_ab wedge R^ab. `internal_covariant` raises the form pair with the Lorentz metric. `legacy_frame_sum` deliberately omits that metric only for comparison. `CS_slice`, `CS_reference`, `relative_slice` and `Q_slice` are per reference volume. Multiply by the appropriate volume to obtain the integral.

Static restriction means H=bD=etaD=0 and b=−V/3. A smooth auxiliary shell sets eta=etaD=0 and b=H, bD=DH. Merely setting an amplitude to zero at one instant does not set its derivative to zero. Nil and Sol reduced trajectories fail the omitted tracefree Einstein equation.

Conversions for spherical scale, torsion inputs, form-Hodge normalization and the frame sum are given in the [convention comparison](legacy_conventions.md).

日本語：`eta,V` は幾何学的torsionの振幅であり、旧入力配列の振幅とは符号が異なる。`H,bD,etaD` は固有時間の量、`qd,ed,Vd` は座標時間微分である。slice量は参照体積あたりの値を保存する。一般の記号式は固定class A族の局所恒等式で、四つの積分値は指定コンパクト商を用いる。完全真空解や一般の異方的解を判定するAPIではない。

## Shape, invariants and endpoint representatives / 追加結果の規約

- The original curvature blocks apply to the common-scale, two-amplitude family. Shape variation starts separately from the LC action after solving the independent connection equation. $`\mathsf B`$ is real symmetric and tracefree; $`B_{\mathrm{end}}`$ is the endpoint generator.
- $`n`$ is the reference matrix, $`m=q^{-1}e^{\mathsf B}ne^{\mathsf B}`$ the physical one. $`X(n)=2n^2-(\mathrm{tr}n)n+4\chi I`$; $`J_2=\mathrm{tr}X^2`$, $`J_3=\mathrm{tr}X^3`$ use unit-scale reference $`n`$ unless stated otherwise.
- EC contractions do not assume pair exchange. $`C_R=R_{ab}{}^{cd}R_{cd}{}^{ef}R_{ef}{}^{ab}`$ sums all ordered indices without bivector factors of one half. The Weyl-square result is restricted to LC.
- Torus linear scale uses $`q=q_0e^{\epsilon\rho}`$; the additive first variation is $`u=q_0\rho`$. Initial shape and shear need not commute. Constant congruence is not a lattice gauge. In `torus_deformations.json`, `a` and `b` are $`a_u`$ and $`b_u`$ of (H-06) under `fixed_shape`, `a` is $`\rho_0`$ of (S-01) and (S-02) under `mixed_shape`, `P` is the Kasner matrix $`\mathsf P`$, and `r` is the coordinate $`\mu=q^{3/2}`$ of (H-07).
- In diagonal evolution $`m_i=n_i e^{2\beta_i}/q`$; $`f,X_i`$ in the ODE use physical $`m`$. The saved ODE writes $`D\beta_i`$ as `s1,s2,s3`. Taylor coefficients use $`X(n)`$ and $`\hat\tau=\tau/q_0`$, saved as `u` in `short_time_evolution.json` and `tau_over_q0` in `short_time_checks.csv`.
- The rotating frame uses integer $`\nu`$ on periods $`2\pi`$ and represents the same static flat metric. The coordinates `x,y,z` in `global_frames.json` are $`x^1,x^2,x^3`$ of the manuscript. It does not identify full dynamical ansatz families or compute spin holonomy.
- After GHY completion, $`P_i=p_i+\partial_iB_{\mathrm{end}}`$ and $`\widetilde p_i=P_i-\partial_iB_{\mathrm{end}}`$. Boundary data must be mapped along with momenta; constant NY adds no bulk dynamics.
- 三つの領域（共通scale・二振幅、torsion消去後のshape、静的frame比較）を区別する。接方向の積分可能性は固定格子の有限斉次sectorの局所的な主張である。Ricci逆問題は一点のスペクトルと構造行列の代数問題であり、局所非等長性・compact quotient全分類は含めない。

The serialized SymPy symbol `theta` is the constant coupling $`\theta_{\mathrm{NY}}`$ of the manuscript.
