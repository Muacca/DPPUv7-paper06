# 付録C. 変分と完全接続方程式

縮約前のscalar密度に(B-02)を代入し、EH加速度項とvector発散を部分積分すると(B-03)となる。体積を掛けたvector項は $`-\mathcal V_{\Sigma}d(q^3V)/\kappa^2`$ であり、HVだけを単独で全微分として捨ててはいけない。

 $`L_0`$ のEuler微分は以下である。ここではFのq微分で $`\eta`$ ,Vを固定する。

式(V-01)。

$$
\begin{aligned}
E_N&=A_{\Sigma}q(F+\dot q^2/N^2),\\
E_{\eta}&=A_{\Sigma}Nq\partial_{\eta}F,\qquad
E_V=A_{\Sigma}Nq\partial_VF,\\
E_q&=A_{\Sigma}\left[N(F+q\partial_qF)+\frac{\dot q^2}{N}
+\frac{2q\ddot q}{N}-\frac{2q\dot q\dot N}{N^2}\right].
\end{aligned}
$$

式(V-01)に補助場解を入れると(H-04)となる。時間変位ξに対して $`\delta q=\xi\dot q,\ \delta\eta=\xi\dot\eta,\ \delta V=\xi\dot V,\ \delta N=\dot{(\xi N)}`$ であり、 $`\dot\xi`$ の項を部分積分すると(H-05)のNoether恒等式が得られる。Nを最初から1に固定すると $`E_N`$ を失うので、gaugeは変分後に選ぶ。

GHYを加えた作用の端点1形式は次のとおりである。これにより(B-06)の固定条件と(B-07)のmomentum shiftが直接得られる。

式(V-02)。

$$
\Theta=\left[p_q\delta q+\delta(B_{\mathrm{EC}}+B_{\mathrm{NY}})\right]_i^f.
$$

完全な接続方程式については、(E-01)のepsilonを逆変換すると $`e^a\wedge T^b-e^b\wedge T^a=0`$ となる。双対frameの内部積を $`\iota_a`$ とし、 $`t=\sum_a\iota_a T^a`$ と置く。aについて内部積を取り足すと、四次元では $`T^b=-e^b\wedge t`$ を得る。さらにbで縮約すると $`t=-3t`$ なので $`t=0`$ 、したがって $`T^b=0`$ である。二振幅への制限を用いないため、除外されたtorsion成分を自由な解として残すことはない。

最後にLCの曲率blockへ $`\eta=V=0`$ を代入すると、 $`\mathrm{Ric}_{00}=-3(DH+H^2)`$ 、 $`\mathrm{Ric}_{ii}=DH+3H^2+(k_j+k_k)/q^2`$ となる。scalarを引けば(E-02)を得る。縮約軌道で $`D^2q=0`$ かつ $`(Dq)^2=-\chi`$ なので、時間成分と空間traceは零となり、 $`\Xi`$ だけが残る。

## 非正則なtorus拘束と端点

補助場消去後のtorusでは $`C=-p_q^2/(4A_{\Sigma}q)`$ で、零集合上で $`dC=0`$ 、Hamiltonベクトル場も零になる。 $`p_q=0`$ は同じ零集合を持つが、 $`q`$ 方向への非零ベクトル場を生成するので同じ拘束gauge系ではない。可逆座標 $`\mu=q^{3/2}`$ 、 $`q>0`$ へ移ると次になる。

式(H-07)。

$$
L=-\frac{4A_{\Sigma}}{9N}\dot\mu^2,\qquad p_\mu=-\frac{8A_{\Sigma}}{9N}\dot\mu,\qquad C=-\frac{9p_\mu^2}{16A_{\Sigma}}.
$$

非正則性は残り、一次の $`\mu=a_\mu+b_\mu\tau`$ にも二次で $`b_\mu=0`$ が必要になる。任意の正の背景lapse $`N_0(t)`$ で $`\tau=\int N_0dt`$ とすれば同じ議論が成立する。端点を固定する $`\xi`$ に対し $`\delta N=d(\xi N_0)/dt`$ なので $`\int\delta Ndt=0`$ 。逆に積分零の変分から $`\xi=N_0^{-1}\int_{t_i}^t\delta Nds`$ を構成できる。固定 $`T`$ のlapse方向はこのgaugeで尽くされるが、 $`q_0`$ のmodulusは除去されない。静的解族の固定時間切片への二形式のpullbackは $`dq_0\wedge dp_q=0`$ であり、この退化を体積modulusのgauge同一視に用いない。

## 端点生成関数と拘束

式(K-01)。

$$
\begin{aligned}
B_{\mathrm{end}}&=-\frac{\mathcal V_{\Sigma}q^3V}{\kappa^2}
 +6\theta_{\mathrm{NY}}\mathcal V_{\Sigma}q^2\eta,\\
P_i&=p_i+\partial_iB_{\mathrm{end}},\qquad
\widetilde p_i=P_i-\partial_iB_{\mathrm{end}},\quad Q^i=(q,N,\eta,V),\\
C'&=-\frac{\widetilde p_q^2}{4A_{\Sigma}q}
 -A_{\Sigma}q(\chi-\eta^2+q^2V^2/9).
\end{aligned}
$$

これはGHY処理後の座標だけに依存する生成関数であり、処理前の速度依存EH端点項を含めない。 $`P_i\,dQ^i=p_i\,dQ^i+dB_{\mathrm{end}}`$ だから二形式は一致する。対称なHessianにより $`\{Q^i,\widetilde p_j\}=\delta^i_j`$ 、 $`\{\widetilde p_i,\widetilde p_j\}=0`$ である。一次拘束は $`\widetilde p_N,\widetilde p_\eta,\widetilde p_V`$ 、二次拘束は $`C^{\prime},\eta,V`$ で、後二つは $`N,q>0`$ で非零係数を除いた形である。

式(K-02)。

$$
\begin{aligned}
\lbrace \zeta_\alpha,\zeta_\beta\rbrace
&=\begin{pmatrix}0&0&-1&0\\
0&0&0&-1\\
1&0&0&0\\
0&1&0&0\end{pmatrix},
&\zeta&=(\widetilde p_\eta,\widetilde p_V,\eta,V),\\
\det\lbrace \zeta_\alpha,\zeta_\beta\rbrace &=1,\qquad
P_\eta=6\theta_{\mathrm{NY}}\mathcal V_{\Sigma}q^2,&
P_V&=-\mathcal V_{\Sigma}q^3/\kappa^2,\\
\Omega_{\mathrm{red}}&=dq\wedge dP_q+dN\wedge dP_N,&
C_{\mathrm{red}}&=-P_q^2/(4A_{\Sigma}q)-A_{\Sigma}q\chi.
\end{aligned}
$$

第二級面では $`\eta=V=0`$ 、 $`P_q=p_q`$ である。二形式のpullbackで $`d\eta=dV=0`$ とすれば上式が直接得られる。 $`\chi=0,P_q=0`$ で $`dC_{\mathrm{red}}=0`$ のままである。固定 $`q,\eta,V`$ は $`\delta B_{\mathrm{end}}=0`$ の十分条件で、lapseの固定は不要である。運動量端点条件を用いる場合は移動式を通じて対応付ける。定数NYの端点代表変更は新しいbulk自由度やtorus拘束の正則化を生じない。

[論文の入口](README_ja.md) / [参考文献](REFERENCES.md)
