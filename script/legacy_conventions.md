# Comparison with earlier conventions

[日本語](legacy_conventions_ja.md) · [Code guide](README.md)

These conversions relate v1 outputs to the geometric conventions of the present paper.

## Spherical scale and torsion map

Let $`\tau`$ denote the old input array and $`T`$ geometric torsion, so $`T=-\tau`$. The former spherical engine used structure magnitude 4 and scale $`\rho`$; the present sphere uses physical radius $`q`$. Their simultaneous map is

Equation (L-01).

$$
\begin{aligned}
\rho&=2q,\qquad \eta_{\tau}=-2\eta,\qquad V_{\tau}=-V,\qquad
p_{\rho}=p_q/2,\\
\mathcal V_{S^3}^{\mathrm{old,corrected}}&=\pi^2/4,
\qquad \mathcal V_{S^3}=2\pi^2,\qquad
(\pi^2/4)\rho^3=2\pi^2q^3.
\end{aligned}
$$

The old spherical value $`\chi=4`$ is correct for its old scale and is compatible with $`\chi=1`$ in the new convention. The uncorrected old code volume $`2\pi^2\rho^3`$ is a separate coefficient error and differs from the corrected volume in (L-01). Keeping the old scales for $`T^3`$, Nil and Sol gives $`\eta_{\tau}=-\eta`$, $`V_{\tau}=-V`$. Their volumes are fixed by the quotients in the text. The same constant factors apply to time derivatives.

## The two duals and the old frame sum

Lower all curvature indices and write $`R_{ab}=\tfrac12 R_{abcd}e^{cd}`$. The former form-Hodge API used the following normalization, related by $`-1/2`$ to the wedge coefficient $`\mathcal P`$ in this paper.

Equation (L-02).

$$
\begin{aligned}
(\ast R_{ab})_{cd}&=\tfrac12\epsilon_{cdef}R_{ab}{}^{ef},\\
P_{\mathrm{form}}^{\mathrm{old}}&=\tfrac14R^{abcd}(\ast R_{ab})_{cd}
=-\tfrac12\mathcal P.
\end{aligned}
$$

The old internal-dual frame sum omitted the factor $`g^{cc}g^{dd}`$ in the $`c\lt d`$ sum. In the same adapted frame, the difference from the covariant quantity is

Equation (L-03).

$$
\mathcal B_{\mathrm{frame}}=3(w\bar y+xz),\qquad
\mathcal B=3(w\bar y-xz),\qquad
\mathcal B_{\mathrm{frame}}-\mathcal B=6xz.
$$

On the static family $`x=z=0`$, so the two expressions agree and (D-03) survives. They differ on general dynamic profiles. For equal $`y_i=y`$ and block values $`(x,y,z,w)=(1,2,3,4)`$, $`\mathcal B=15`$. A 01-boost with $`\cosh\beta=5/4`$ and $`\sinh\beta=3/4`$ keeps it at 15 but changes the frame sum from 33 to $`231/4`$. The old sum is not used as a Lorentz scalar.
