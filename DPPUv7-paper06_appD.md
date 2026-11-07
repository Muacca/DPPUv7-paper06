# 付録D. 斉次shape変分と局所真空発展

## 一般行列の作用と五方向変分

独立接続方程式によりtorsionを消去し、GHY完成したLC作用を使う。固定class-A coframeで $`e=q e^{\mathsf B}\sigma`$ とすると、空間曲率は $`f(n_{\mathsf B})/q^2`$ 、外部曲率のtracefree部分は $`\mathsf S`$ である。ここで $`D=N^{-1}d/dt`$ である。

式(S-04)。

$$
\begin{aligned}
n_{\mathsf B}&=e^{\mathsf B}ne^{\mathsf B},\qquad
\mathsf S=\tfrac12e^{-\mathsf B}D(e^{2\mathsf B})e^{-\mathsf B},\\
L&=\frac{\mathcal V_{\Sigma}Nq^3}{2\kappa^2}
 \left[-6H^2+\mathrm{tr}(\mathsf S^2)+\frac{f(n_{\mathsf B})}{q^2}\right].
\end{aligned}
$$

$`f`$ の式は対角 $`n`$ のKoszul接続から得た空間Ricciのtraceであり、実対称行列の定数回転で一般 $`n`$ へ成立する。 $`\mathsf B,D\mathsf B`$ をともに一次の量と数え、非可換積の順序を保って指数行列を展開する。 $`\mathsf S=D\mathsf B+O(3)`$ で、二次のcommutator項は両側の $`e^{-\mathsf B}`$ により相殺する。ポテンシャルは次になる。

式(S-05)。

$$
\begin{aligned}
f(n_{\mathsf B})&=6\chi-2\mathrm{tr}(X\mathsf B)+U_2+O(3),\\
U_2&=2[\mathrm{tr}(\mathsf Bn)]^2
 +2(\mathrm{tr}n)\mathrm{tr}(\mathsf B^2n)\\
&\hspace{1em}-4\mathrm{tr}(\mathsf Bn\mathsf Bn)
 -4\mathrm{tr}(\mathsf B^2n^2).
\end{aligned}
$$

実際、 $`n_{\mathsf B}=n+\mathsf Bn+n\mathsf B+\mathsf B^2n/2+n\mathsf B^2/2+\mathsf Bn\mathsf B+O(3)`$ を $`f`$ に代入すれば各項が得られる。 $`\mathrm{tr}\mathsf B=0`$ より一次項は式(E-04)となる。五方向は $`\mathsf B_{11}=\mathsf b_1,\mathsf B_{22}=\mathsf b_2,\mathsf B_{33}=-\mathsf b_1-\mathsf b_2`$ と独立な三非対角成分で尽くされる。その速度Gram形式は $`\dot{\mathsf b}_1^2+\dot{\mathsf b}_2^2+(\dot{\mathsf b}_1+\dot{\mathsf b}_2)^2+2(\dot{\mathsf b}_3^2+\dot{\mathsf b}_4^2+\dot{\mathsf b}_5^2)`$ で正定値である。

$`X=0`$ を分類するため $`n`$ を対角化する。差 $`X_i-X_j=(n_i-n_j)[2(n_i+n_j)-s_1]`$ は、三つの相異なる固有値を許さない。二つが等しい場合 $`n=(a,a,b)`$ とおけば差は $`b(b-a)`$ に比例し、 $`b=a`$ または $`b=0`$ となる。各枝を式(E-04)へ戻すと実際に $`X=0`$ である。これで式(E-05)の必要十分性が証明される。

共通scale背景まわりのmomentum方程式も独立に評価すると、 $`\delta G_{01}=(n_2-n_3)D\mathsf B_{23}/q`$ と巡回置換を得る。したがって一般の非対角shape速度を無条件には許さない。torusでは $`n=0`$ なので全て零、対角class-A sectorでは非対角速度が零なので満たされる。

## Torusのlapse拘束と厳密解族

$`n=0`$ ではポテンシャルが恒等的に零である。 $`q=q_0e^{\epsilon\rho}`$ 、 $`\mathsf B=\epsilon\mathsf B_1`$ 、 $`N=N_0+\epsilon\delta N`$ を式(S-04)へ代入する。二次項は式(S-01)、三次項と元のlapse方程式は次を与える。

式(S-06)。

$$
\begin{aligned}
L^{[3]}&=\frac{A_{\Sigma}q_0^3}{6N_0}\left(3\rho-\frac{\delta N}{N_0}\right)
 [\mathrm{tr}(\dot{\mathsf B}_1^2)-6\dot\rho^2],\\
E_N^{[2]}&=\frac{A_{\Sigma}q_0^3}{6N_0^2}[6\dot\rho^2-\mathrm{tr}(\dot{\mathsf B}_1^2)]
 =\partial_{\delta N}L^{[3]}.
\end{aligned}
$$

二次作用にlapse一次変分が現れないことは拘束の消失を意味しない。元のlapse変分を先に取るか、三次作用を $`\delta N`$ で変分して同じ二次拘束を得る。二次Euler式は固有時間で $`D^2\rho=D^2\mathsf B_1=0`$ なので、式(S-01)のaffine解と錐が必要である。

十分性は式(S-02)の厳密族で示す。定数congruence前の計量で混合外部曲率を $`K=(3\epsilon h/\psi)\mathsf P`$ とすると、空間Ricciは零、 $`DK+(\mathrm{tr}K)K=0`$ である。時間Ricciは次の通りで、錐上では零になる。

式(S-07)。

$$
\mathrm{Ric}_{00}=\frac{9\epsilon^2h^2}{\psi^2}[1-\mathrm{tr}(\mathsf P^2)]=-\frac{\epsilon^2}{\psi^2}[\mathrm{tr}(\Sigma^2)-6h^2].
$$

定数congruenceは局所Ricci零を保ち、各時刻の定数係数正定値計量は同じ格子上に定義できる。これは格子変換をgaugeとする主張ではない。行列指数の一次展開は $`g_\epsilon=q_0^2[I+2\epsilon(\rho_0I+\mathsf B_0+h\tau I+\Sigma\tau)]+O(\epsilon^2)`$ となり、初期shapeと $`\Sigma`$ が非可換でも指定接方向を得る。 $`\det g_\epsilon`$ から式(S-02)の $`q_\epsilon`$ が従う。 $`h=0`$ 枝は正定値性から $`\Sigma=0`$ となり、静的族で実現する。

任意の一次lapse変分は正の $`N_\epsilon=N_0+\epsilon\delta N+O(\epsilon^2)`$ と $`\tau_\epsilon=\int N_\epsilon dt`$ を選び、厳密族へ合成して実現できる。背景が静的なのでこの合成は一次の空間計量接方向を変更しない。局所・初期値の比較と、両端計量を固定する変分問題は区別する。後者では各affine係数が両端で零となり、全ての一次scale・shape係数が零になる。

## 対角sectorの閉性と短時間係数

対角 $`n`$ と $`e^i=qe^{\beta_i}\sigma^i`$ 、 $`\sum_i\beta_i=0`$ を用いる。 $`m_i=n_ie^{2\beta_i}/q`$ とし、以後の $`f,X_i`$ は物理行列 $`m`$ で評価する。Koszul接続の対角空間Ricci固有値は $`2m_i^2-(\sum m)m_i+f`$ で、そのtracefree部分は $`X_i=2m_i^2-(\sum m)m_i+2f/3`$ である。外部曲率は対角で、momentumと非対角Einstein方程式は零になるので、このsectorは閉じている。

式(S-08)。

$$
\begin{aligned}
Dq&=qH,& DH&=-3H^2-f/3,\\
D^2\beta_i&=-3HD\beta_i-X_i(m),& Dm_i&=m_i(2D\beta_i-H),\\
\mathcal C&=f+6H^2-\sum_i(D\beta_i)^2=0,&D\mathcal C&=-6H\mathcal C.
\end{aligned}
$$

拘束伝播は $`Df=-2Hf-2\sum_iX_iD\beta_i`$ を使って式(S-08)を微分すれば直接確認できる。また拘束上で $`D^2q=-q\sum_i(D\beta_i)^2/3`$ である。初期 $`q=q_0,\beta=D\beta=0,H_0^2=-\chi/q_0^2`$ では $`D^2\beta_i=-X_i(n)/q_0^2`$ 、 $`DX_i(m)=-2H_0X_i(n)/q_0^2`$ だから $`D^3\beta_i=5H_0X_i(n)/q_0^2`$ となる。 $`D^2q=D^3q=0`$ と次の四次係数が得られ、式(S-03)を証明する。

式(S-09)。

$$
D^4q\big|_0=-\frac{2q_0}{3}\sum_i(D^2\beta_i|_0)^2=-\frac{2J_2(n)}{3q_0^3}.
$$

補足資料は $`q_0=1`$ 、Nil・Sol、 $`H_0`$ の両符号で $`\hat\tau=0.01,0.02,0.04`$ の計12点をDOP853で対照する。指定値は rtol $`2\times10^{-13}`$ 、atol $`2\times10^{-15}`$ 、最大step $`0.002`$ である。拘束残差とTaylorとの差を保存し、検査点の対照と局所Taylor次数を区間誤差保証と区別する。

[論文の入口](README_ja.md) / [参考文献](REFERENCES.md)
