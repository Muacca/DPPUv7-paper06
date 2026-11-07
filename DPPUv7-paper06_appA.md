# 付録A. 構造定数と成分規約

## JacobiとBianchi class

構造定数を $`d\sigma^i=\tfrac12C^i{}_{jk}\sigma^j\wedge\sigma^k`$ で定め、五係数で表示する。scaleを外した非零成分は以下であり、下の二添字は反対称とする。Jacobi残差とtraceにより、この表示がLie代数を定める条件とunimodular条件が得られる。

式(A-01)。

$$
\begin{aligned}
C^1{}_{23}&=a,\quad C^2{}_{31}=b,\quad C^3{}_{12}=c,\quad
C^2{}_{12}=u,\quad C^3{}_{13}=v,\\
J_{123}&=(-a(u+v),0,0),\qquad C^j{}_{ij}=(u+v,0,0).
\end{aligned}
$$

したがってLie familyには $`a(u+v)=0`$ が必要であり、unimodularな主線では $`u+v=0`$ であり、これは式(G-01)で $`n=\mathrm{diag}(a,b,c)`$ 、 $`n_{23}=n_{32}=u`$ とした場合である。Jacobi条件外の多項式恒等式をLie群の曲率命題と読まない。Solの $`u=1,v=-1`$ はclass Aであり、23平面の向きを保つπ/4回転で構造行列を対角化できる。本文の一般式はこのclass A主線に対するものである。

## 完全NYの成分規約

完全NYの成分係数は、向き付けた交代記号[0123]＝＋1で書くとTT側が1/4、曲率側が1/2となる。この交代記号を、添字を上げたLorentz epsilonと混同しない。

式(A-02)。

$$
\mathcal N=\tfrac14[abcd]T^e{}_{ab}T_{ecd}
-\tfrac12[abcd]R_{abcd}.
$$

## 固定torus上の周期回転coframe

$`x^1,x^2,x^3`$ の周期を $`2\pi`$ とし、 $`\nu\in\mathbb Z`$ とする。次の $`O(x^3)\in SO(3)`$ は大域的であり、 $`\sigma=O(dx^1,dx^2,dx^3)^T`$ は同じflat計量の正規直交coframeである。

式(F-01)。

$$
\begin{aligned}
O(x^3)&=\begin{pmatrix}\cos\nu x^3&\sin\nu x^3&0\\
-\sin\nu x^3&\cos\nu x^3&0\\
0&0&1\end{pmatrix},\\
d\sigma^1&=-\nu\sigma^2\wedge\sigma^3,\quad
d\sigma^2=-\nu\sigma^3\wedge\sigma^1,\quad d\sigma^3=0,\\
n&=\mathrm{diag}(-\nu,-\nu,0),\qquad
\omega=-dO\thinspace O^T,\quad \omega^1{}_2=-\nu dx^3.
\end{aligned}
$$

直接計算で $`d\sigma+\omega\wedge\sigma=0`$ 、 $`d\omega+\omega\wedge\omega=0`$ となる。周期整数条件から $`O(x^3+2\pi)=O(x^3)`$ であり、vector平行移動は各閉周期で恒等になる。 $`\nu\ne0`$ のとき、定数 $`\alpha_i`$ に対して $`d(\alpha_i\sigma^i)=0`$ は $`\alpha_1=\alpha_2=0`$ を与える。一方 $`dx^1,dx^2,dx^3`$ はこのcoframeでは可変係数を持つ大域閉形式で、その三周期に対するperiod行列は $`2\pi I`$ である。

$`b_1=3`$ の上限もflat torusのFourier分解から得られる。閉1形式の非零波数成分はその波数ベクトルに比例するためexactであり、cohomologyには零波数の三成分だけが残る。したがって定数係数coframe部分空間の閉形式数はBetti数ではない。これは指定した静的計量とvector接続の比較であり、spin構造や両frameの動的ansatz全体の同値性を述べない。

[論文の入口](README_ja.md) / [参考文献](REFERENCES.md)
