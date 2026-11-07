# 3. Bulk作用と符号付き境界完成

同じ接続の曲率を用い、完全NYを外微分として作用に入れる。 $`\theta_{\mathrm{NY}}`$ は定数、 $`\kappa^2`$ は正とする。境界は時間区間の両端にある空間的sliceであり、空間商は閉じている。

式(B-01)。

$$
\begin{aligned}
S&=\int_M e^{0123}\frac{R_{\mathrm{EC}}}{2\kappa^2}
+\theta_{\mathrm{NY}}\int_M dQ+S_{\mathrm{GHY}},\\
R^a{}_{b}&=d\omega^a{}_{b}+\omega^a{}_{c}\wedge\omega^c{}_{b},\qquad Q=e^a\wedge T_a,\\
dQ&=T^a\wedge T_a-e^a\wedge e^b\wedge R_{ab}.
\end{aligned}
$$

接続を時間依存のまま微分すると、scalar曲率は次のように分解する。vector振幅の微分項とHとの積は、体積を掛けたときに一つの発散をなす。

式(B-02)。

$$
\begin{aligned}
R_{\mathrm{LC}}&=6(DH+2H^2+\chi/q^2),\\
R_{\mathrm{EC}}-R_{\mathrm{LC}}&=-6\eta^2/q^2+2V^2/3-2DV-6HV.
\end{aligned}
$$

式(B-03)。

$$
\begin{aligned}
L_{\mathrm{raw}}&=L_0+\frac d{dt}(B_{\mathrm{EH}}+B_{\mathrm{EC}}+B_{\mathrm{NY}}),\\
L_0&=A_{\Sigma}q\left(NF-\frac{\dot q^2}{N}\right),\qquad
A_{\Sigma}=\frac{3\mathcal V_{\Sigma}}{\kappa^2},\\
F&=\chi-\eta^2+\frac{q^2V^2}{9}.
\end{aligned}
$$

式(B-04)。

$$
B_{\mathrm{EH}}=\frac{3\mathcal V_{\Sigma}q^2\dot q}{\kappa^2N},\qquad
B_{\mathrm{EC}}=-\frac{\mathcal V_{\Sigma}q^3V}{\kappa^2},\qquad
B_{\mathrm{NY}}=6\theta_{\mathrm{NY}}\mathcal V_{\Sigma}q^2\eta.
$$

式(B-03)は四幾何に共通である。構造行列はFに $`\chi`$ を通じて入り、商の体積は $`A_{\Sigma}`$ へ入る。完全NYの寄与は $`B_{\mathrm{NY}}`$ に集約されるため、bulkの補助場方程式と縮約軌道は定数NY結合によらない。

外向き単位法線は最終端で正、初期端で負の時間方向を向く。LC接続で定義した外部曲率とsignature因子を含むGHY項を用いる[GH1977,HH1996]。

式(B-05)。

$$
\begin{aligned}
\hat n_f&=+N^{-1}\partial_t,\quad \hat n_i=-N^{-1}\partial_t,\quad
\epsilon_n=\hat n^a\hat n_a=-1,\quad K_f=3H_f,\quad K_i=-3H_i,\\
S_{\mathrm{GHY}}&=\frac1{\kappa^2}\sum_{s=i,f}\epsilon_n\int_{\Sigma_s}\sqrt h K_s
=-[B_{\mathrm{EH}}]_i^f.
\end{aligned}
$$

式(B-06)。

$$
S=\int_{t_i}^{t_f}L_0dt+[B_{\mathrm{EC}}+B_{\mathrm{NY}}]_i^f,\qquad
\delta q=\delta\eta=\delta V=0\quad\text{at both endpoints}.
$$

この固定端点条件は完成作用の変分に十分であり、端点のlapseやscale速度を固定する必要はない。 $`L_0`$ をcanonical代表に選ぶと、残る全微分B＝ $`B_{\mathrm{EC}}`$ ＋ $`B_{\mathrm{NY}}`$ を戻した代表の運動量は次の分だけ変わる。bulk方程式は変わらない。

式(B-07)。

$$
\begin{aligned}
\Delta p_q&=-\frac{3\mathcal V_{\Sigma}q^2V}{\kappa^2}
+12\theta_{\mathrm{NY}}\mathcal V_{\Sigma}q\eta,\\
\Delta p_{\eta}&=6\theta_{\mathrm{NY}}\mathcal V_{\Sigma}q^2,\qquad
\Delta p_V=-\frac{\mathcal V_{\Sigma}q^3}{\kappa^2},\qquad \Delta p_N=0.
\end{aligned}
$$

端点で非零torsionを指定することと、その端点値を持つ滑らかな停留解が存在することは別である。第4節の補助場方程式を区間全体に課すと、連続性により端点でも $`\eta=V=0`$ となる。境界sourceや不連続な層を追加する問題はここでは扱わない。

## 端点代表の正準対応

GHY処理後の端点関数を $`B_{\mathrm{end}}=B_{\mathrm{EC}}+B_{\mathrm{NY}}`$ と書く。 $`L_0+dB_{\mathrm{end}}/dt`$ は $`P_i=p_i+\partial_iB_{\mathrm{end}}`$ を与え、補助変数の運動量も同時に移動する。 $`\widetilde p_i=P_i-\partial_iB_{\mathrm{end}}`$ を用いると拘束は元と同じ形になる。付録Cで二形式と第二級拘束を示す。正準同値性は対応する端点データについての主張であり、旧運動量固定と新運動量固定をそのまま同じ境界条件とはしない。

[論文の入口](README_ja.md) / [参考文献](REFERENCES.md)
