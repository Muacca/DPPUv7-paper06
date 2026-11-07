# 付録B. 曲率恒等式の導出

対角化した構造方程式では $`de^i=He^{0i}+(n_i/q)e^{jk}`$ となる。式(G-04)の空間係数は次の関係を満たすので、第一Cartan方程式の空間torsionは2 $`\eta`$ /qとなる。

式(R-01)。

$$
a_j+a_k=\frac{n_i-2\eta}{q},\qquad
\frac{n_i}{q}-a_j-a_k=\frac{2\eta}{q},\qquad H-b=V/3.
$$

第二Cartan方程式を展開すると、boost曲率の時間係数は $`Db+Hb`$ 、空間係数は $`b(n_i/q-a_j-a_k)`$ である。回転曲率の時間係数と空間係数は以下となる。 $`\eta`$ に一次の空間項が相殺する点が共通公式の鍵である。

式(R-02)。

$$
\begin{aligned}
Da_i+Ha_i&=-D\eta/q,\\
b^2+\frac{n_i a_i}{q}-a_ja_k
&=b^2+\frac{n_i\ell_i-\ell_j\ell_k-\eta^2}{q^2},\\
\ell_j+\ell_k&=n_i.
\end{aligned}
$$

これにより(D-01)が得られる。空間scalarとRicciはtorsion、b、Dbを零にした曲率から導く。 $`k_i`$ を足すと(G-02)となり、四代表点の $`\chi`$ とtracefree Ricciが得られる。任意の対称nは一定SO(3)回転でこのframeへ移せるので、スカラー恒等式は対角成分だけに依存する表示の人工物ではない。

全内部添字を足したPontryagin外積では、三つのboost対が $`-12xw`$ 、回転対が $`4z\sum_i y_i`$ を与える。内部双対縮約では、空間形式に重み＋1、時間を含む形式に重み−1が掛かるため、 $`w\sum_i y_i-3xz`$ となる。したがって(D-02)が成立する。

静的制限 $`H=Db=D\eta=0`$ 、 $`b=-V/3`$ を取ると、 $`\mathcal B=3w\bar y`$ となる。この多項式のVηに一次の係数を抽出すると $`C=-9\chi`$ を得る。

NYは直接外積でも確認できる。両側のηV項が相殺し、Qの時間微分だけが残る。

式(R-03)。

$$
\begin{aligned}
T^a\wedge T_a&=\frac{4V\eta}{q}e^{0123},\\
e^a\wedge e^b\wedge R_{ab}
&=\left(\frac{4V\eta}{q}-\frac{6(D\eta+2H\eta)}q\right)e^{0123}.
\end{aligned}
$$

CS形式を接続から展開すると、 $`\eta=0`$ の定数部分は次の多項式となる。基本対称式へ整理すると(D-05)の $`\Phi_0`$ である。

式(R-04)。

$$
\Phi_0=\sum_i n_i^3-\sum_{i\ne j}n_i^2n_j+4n_1n_2n_3
=s_1^3-4s_1s_2+10s_3.
$$

残る $`\phi`$ を微分し、 $`Dq=Hq`$ と $`Db=x-Hb`$ を用いると、 $`-12(\chi-\eta^2)D\eta-12q^2b^2D\eta-24q^2b\eta x=q^3\mathcal P`$ となる。これは局所密度と端点primitiveの独立な照合でもある。

一般の相対形式(D-06)は接続の経路 $`\omega_\gamma=\bar\omega+\gamma K`$ から導ける。 $`R_\gamma=\bar R+\gamma\bar D K+\gamma^2K\wedge K`$ を用い、Bianchi恒等式とtraceの巡回性で次式を得る。 $`\gamma`$ を0から1へ積分すると(D-06)となる。この証明は特定の四代表点に依存しない。

式(R-05)。

$$
\partial_\gamma[-\mathrm{tr}(R_\gamma\wedge R_\gamma)]
=-2d\thinspace \mathrm{tr}(K\wedge R_\gamma).
$$

## 全添字縮約とLC Weyl二乗

$`X`$ の固有値を $`X_i`$ とすると $`y_i=\bar y-X_i/q^2`$ 、 $`\sum_iX_i=0`$ である。従って $`\sum y_i=3\bar y`$ 、 $`\sum y_i^2=3\bar y^2+J_2/q^4`$ 、 $`\sum y_i^3=3\bar y^3+3\bar yJ_2/q^4-J_3/q^6`$ となる。曲率blockの添字を明示して和を取ると次を得る。

式(I-06)。

$$
\begin{aligned}
I_1&=4\sum_i(x^2+y_i^2-w^2-z^2),\\
I_2&=4\sum_i(x^2+y_i^2+2wz),\\
C_R&=8\sum_i[x^3+y_i^3+3wz(x+y_i)],\\
\mathrm{Ric}^a{}_b&=\mathrm{diag}
 (3x,x+2\bar y+X_1/q^2,x+2\bar y+X_2/q^2,x+2\bar y+X_3/q^2).
\end{aligned}
$$

前二式の負号と交差項の違いはLorentz計量と添字対交換の違いによる。最後の対角行列を二乗・三乗してtraceを取ると式(I-02)、(I-03)のRicci縮約になる。LCでは $`w=z=0`$ 、 $`R=6(x+\bar y)`$ である。四次元のtrace除去式 $`C^2=I_1-2I_{\mathrm{Ric}}+R^2/3`$ に代入すると、 $`x,\bar y`$ の項がすべて消え、式(I-04)を得る。これらの計算は選定した縮約についての完結した導出である。

## Ricciスペクトルと構造行列の逆問題

実対称tracefree $`X`$ についてNewton恒等式は $`\sum_{i\lt j}X_iX_j=-J_2/2`$ 、 $`\det X=J_3/3`$ を与える。よって式(I-05)の特性多項式が得られる。三つの実根を持つため判別式 $`J_2^3/2-3J_3^2`$ は非負である。逆に $`J_2\ge0`$ とこの判別式条件は三実根を与える。 $`J_2=0`$ は $`X=0`$ 、判別式零は重複根を表すが、後述のRicci零固有値とは別の条件である。

単位scaleで $`r_i=X_i+2\chi`$ を空間Ricci固有値とする。 $`s_1=\mathrm{tr}n`$ 、 $`\ell_i=s_1/2-n_i`$ とおくと、 $`n_i=\ell_j+\ell_k`$ と $`r_i=2\ell_j\ell_k`$ が成立する。これは $`k_j+k_k`$ の式へ代入して直接得られる。実際のscaleではRicci固有値は $`r_i/q^2`$ である。

式(I-07)。

$$
r_i\ne0\ (\forall i):\qquad r_1r_2r_3>0,\qquad \ell_i^2=\frac{r_jr_k}{2r_i},\qquad n_i=\ell_j+\ell_k.
$$

非退化枝では積の正値性が必要であり、式(I-07)の正の平方根の大きさと $`r_i=2\ell_j\ell_k`$ の符号条件から十分性も従う。例えば $`\ell_1`$ の符号を選ぶと $`\ell_2=r_3/(2\ell_1)`$ 、 $`\ell_3=r_2/(2\ell_1)`$ が定まり、残る式は自動的に満たされる。したがって全体符号と固有空間の直交基底選択を除き一意である。実対称行列の置換・共役は $`SO(3)`$ 共役として実現できる。

二つの $`r_i`$ が非零なら三つの $`\ell_i`$ がすべて非零となり、残る $`r_i`$ も非零なので「ちょうど二つ」は不可能である。 $`r_1\ne0,r_2=r_3=0`$ なら $`\ell_1=0`$ 、 $`\ell_2\ell_3=r_1/2`$ で、積だけが固定される。例えば $`\ell=(0,t,t^{-1})`$ 、 $`t\ne0`$ は $`n=(t+t^{-1},t^{-1},t)`$ と同じ $`r=(2,0,0)`$ を与える。すべての $`r_i`$ が零なら、非零 $`\ell_i`$ は高々一つであり、 $`n\sim(a,a,0)`$ となる。これは構造行列の代数的逆問題であり、局所等長類やcompact quotientの分類ではない。

[論文の入口](README_ja.md) / [参考文献](REFERENCES.md)
