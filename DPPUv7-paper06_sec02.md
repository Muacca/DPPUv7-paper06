# 2. 幾何と変分領域

向き付けた正規直交coframeを用いる。空間構造行列nは実対称で時間によらず、空間frameの一定SO(3)回転で対角化できる。空間添字は1,2,3、四次元添字は0,1,2,3であり、構造方程式の添字は総和する。

式(G-01)。

$$
\begin{aligned}
g_{ab}&=\mathrm{diag}(-1,1,1,1),\qquad \epsilon_{0123}=+1,\quad \epsilon^{0123}=-1,\\
e^0&=Ndt,\qquad e^i=q\sigma^i,\qquad
d\sigma^i=\tfrac12 n_{ij}\epsilon_{jkl}\sigma^k\wedge\sigma^l,\\
D&=N^{-1}d/dt,\qquad H=Dq/q,\qquad
\mathcal V_{\Sigma}=\int_{\Sigma}\sigma^{123},\quad q,N>0.
\end{aligned}
$$

このclass A表示はJacobi恒等式を満たす。対角frameで固有値を $`n_i`$ とし、参照空間のsectional曲率係数を $`k_i`$ と書く。空間Ricciのi方向は残る二つの $`k_i`$ の和である。 $`\chi`$ の値はこの曲率から定義する。

式(G-02)。

$$
\begin{aligned}
\ell_i&=\tfrac12(n_1+n_2+n_3)-n_i,\qquad
k_i=n_i\ell_i-\ell_j\ell_k,\\
{}^{(3)}R&=\frac{2(k_1+k_2+k_3)}{q^2}=\frac{6\chi}{q^2},\\
6\chi&=\tfrac12(\mathrm{tr}\thinspace n)^2-\mathrm{tr}(n^2),\qquad
{}^{(3)}\mathrm{Ric}_{ii}=\frac{k_j+k_k}{q^2}.
\end{aligned}
$$

式(G-02)のi,j,kは循環する三組を表し、追加の総和は取らない。四つの参照幾何を以下に固定する。球面のqは物理半径であり、Paper05の規約と同じである。

| 幾何 | nの固有値 |  $`\chi`$  | 参照体積 |
| --- | --- | --- | --- |
|  $`T^3`$  | $`(0,0,0)`$ | 0 | $`(2\pi)^3`$ |
| Nil | $`(0,0,-1)`$ | $`-1/12`$ | $`(2\pi)^4`$ |
| Sol | $`(0,1,-1)`$ | $`-1/3`$ | $`\log((3+\sqrt5)/2)`$ |
|  $`S^3`$  | $`(2,2,2)`$ | 1 | $`2\pi^2`$ |

torusの三つの角座標の周期は2πとする。空間座標を $`x^1,x^2,x^3`$ と書く。Nilではcoframeを $`dx^1,dx^2,dx^3-x^1dx^2`$ とし、格子生成元を $`(x^1,x^2,x^3)\mapsto(x^1+2\pi,x^2,x^3+2\pi x^2)`$ 、 $`(x^1,x^2,x^3)\mapsto(x^1,x^2+2\pi,x^3)`$ 、 $`(x^1,x^2,x^3)\mapsto(x^1,x^2,x^3+4\pi^2)`$ とする。いずれの同定もcoframeを保つ。

Solではcoframeを $`dx^1,e^{x^1}dx^2,e^{-x^1}dx^3`$ とする。A＝((2,1),(1,1))の固有値λ＝(3+√5)/2とλ⁻¹を使い、 $`x^1`$ を $`\log\lambda`$ だけ進めるとfiber座標をdiag(λ⁻¹,λ)で同定する。fiber格子はAの固有座標へ移した面積1の格子である。coframeはこのmapping torusへ降りる。構造行列は対角化前には23成分と32成分が1であり、class AのVI₀である[paper4]。

幾何学的torsionを第一Cartan方程式で定義し、二つの実振幅を次式で固定する。ここでも空間外積は循環する各組について記す。

式(G-03)。

$$
T^a=de^a+\omega^a{}_{b}\wedge e^b,\qquad
T^0=0,\qquad T^i=\frac V3e^{0i}+\frac{2\eta}{q}e^{jk}.
$$

式(G-04)。

$$
\begin{aligned}
a_i&=(\ell_i-\eta)/q,\qquad b=H-V/3,\\
\omega^0{}_{i}&=\omega^i{}_{0}=b e^i,\qquad
\omega^j{}_{k}=a_i e^i,\qquad \omega^k{}_{j}=-a_i e^i.
\end{aligned}
$$

式(G-04)を(G-03)へ代入すると指定torsionが再現される。時間profileは十分滑らかとし、空間計量の形状と商の格子は固定する。一般の構造行列に対する式は局所恒等式であり、具体的な空間積分には上記の固定コンパクト商を用いる。scaleとtorsionの変分を先に制限した方程式と、独立接続・coframeの完全方程式との関係は第5節で検査する。

## 変分領域の区別

共通scale・二振幅族では固定した実対称構造行列 $`n`$ と $`q,N,\eta,V`$ を用いる。第6節の時間依存曲率公式の適用域はこの族である。shape解析では独立接続方程式を先に解き、LC作用で $`e^i=q(e^{\mathsf B})^i{}_j\sigma^j`$ 、 $`\mathsf B^{T}=\mathsf B`$ 、 $`\mathrm{tr}\mathsf B=0`$ とする。 $`\mathsf B`$ は五成分の時間依存行列であり、端点生成関数 $`B_{\mathrm{end}}`$ とは異なる。物理構造行列は $`m=q^{-1}e^{\mathsf B}ne^{\mathsf B}`$ である。

固定した大域coframeに定数の正定値係数を与えた計量は同じ商上に定義できる。定数congruenceを格子保存座標変換と同一視しない。第7節・付録Aの回転coframeは同じ静的flat torus計量の表示比較であり、それぞれのframeに対する動的な定数係数shape族の同一性は仮定しない。

[論文の入口](README_ja.md) / [参考文献](REFERENCES.md)
