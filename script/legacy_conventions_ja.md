# 旧規約との比較

[English](legacy_conventions.md) · [計算案内](README_ja.md)

v1の出力と本稿の幾何学的規約を比較するため、半径・torsion振幅・双対縮約の対応を示す。

## 球面とtorsionの換算

旧実装の入力配列を $`\tau`$ 、幾何学的torsionをTとすると $`T=-\tau`$ である。旧球面の構造係数の大きさ4とscaleρに対し、新しい球面は物理半径qを使う。両者の同時換算は次式となる。

式(L-01)。

$$
\begin{aligned}
\rho&=2q,\qquad \eta_{\tau}=-2\eta,\qquad V_{\tau}=-V,\qquad
p_{\rho}=p_q/2,\\
\mathcal V_{S^3}^{\mathrm{old,corrected}}&=\pi^2/4,
\qquad \mathcal V_{S^3}=2\pi^2,\qquad
(\pi^2/4)\rho^3=2\pi^2q^3.
\end{aligned}
$$

旧球面の $`\chi=4`$ は旧scaleの定義に対して正しく、新規約の $`\chi=1`$ と矛盾しない。旧コードの未訂正体積 $`2\pi^2\rho^3`$ は別の係数誤りであり、式(L-01)の訂正体積と区別する。 $`T^3`$ 、Nil、Solの旧scaleを維持するとtorsion振幅の対応は $`\eta_{\tau}=-\eta`$ 、 $`V_{\tau}=-V`$ となる。体積は本文の固定商に合わせる。時間微分にも同じ定数因子を適用する。

## 二つの双対と旧frame和

曲率成分をすべて下げて $`R_{ab}=\tfrac12 R_{abcd}e^{cd}`$ とする。旧form-Hodge APIの規格化は次式であり、本稿の外積係数Pとは−1/2の関係にある。

式(L-02)。

$$
\begin{aligned}
(\ast R_{ab})_{cd}&=\tfrac12\epsilon_{cdef}R_{ab}{}^{ef},\\
P_{\mathrm{form}}^{\mathrm{old}}&=\tfrac14R^{abcd}(\ast R_{ab})_{cd}
=-\tfrac12\mathcal P.
\end{aligned}
$$

内部双対の旧frame和では、c＜dの和に必要な $`g^{cc}g^{dd}`$ が欠けていた。同じ適応frameにおける旧和と共変な量の差は次のとおりである。

式(L-03)。

$$
\mathcal B_{\mathrm{frame}}=3(w\bar y+xz),\qquad
\mathcal B=3(w\bar y-xz),\qquad
\mathcal B_{\mathrm{frame}}-\mathcal B=6xz.
$$

静的族では $`x=z=0`$ なので両者は一致し、式(D-03)を保持できる。一般の動的族では一致しない。例えば等しい三つのyとblock値 $`(x,y,z,w)=(1,2,3,4)`$ を取ると $`\mathcal B=15`$ である。 $`\cosh\beta=5/4`$ 、 $`\sinh\beta=3/4`$ の01-boostでBは15に保たれるが、frame和は33から231/4へ変化する。旧和をLorentz scalarとして用いない。
