## 2. Lorentzian-native setup and topology conventions

### 2.1 Reduced Lorentzian frame

本稿では Lorentzian 符号 $(-,+,+,+)$ の orthonormal frame を用いる。coframe は

$$
e^0 = N(t)\,dt, \qquad e^i = q(t)\,\sigma^i_\Sigma
$$

と書く。 $N(t)$ は lapse function、 $q(t)$ は topology $\Sigma$ に共通の isotropic scale、 $\sigma^i_\Sigma$ は $\Sigma$ の左不変 coframe である。torsion 変数は axial component $\eta(t)$ および vector-trace component $V(t)$ を用い、それぞれ AX および VT の pure branch、混合した MX branch を与える。

空間 3-curvature の正規化を

$$
{}^{(3)}R = \frac{6\chi}{q^2}
$$

で定義し、 $\chi$ を本稿全体の topology scalar として使う。 $\chi$ は EH（pure Einstein-Hilbert）sector からのみ導かれ、torsion 変数には依存しない。

### 2.2 Topology structure constants and $\chi$ values

四つの topology の構造定数は、Bianchi 型 homogeneous space の class-A / class-B 分類 [11] および本 DPPU 系列の三 topology mode dictionary [4] の枠組みに沿って、five-parameter LZ-native isotropic family

$$
C^1{}_{23}=\frac{a}{q},\quad C^2{}_{31}=\frac{b}{q},\quad C^3{}_{12}=\frac{c}{q},\quad C^2{}_{12}=\frac{u}{q},\quad C^3{}_{13}=\frac{v}{q}
$$

の特殊化として与えられる（下添え字は反対称）。このとき

$$
{}^{(3)}R\,q^2=
-\tfrac{1}{2}\!\left(a^2-2ab-2ac+b^2-2bc+c^2+4u^2+4uv+4v^2\right)
=:-\tfrac{1}{2}\,\Delta,
$$

である。定義 ${}^{(3)}R=6\chi/q^2$ と比較すると

$$
6\chi=-\frac{\Delta}{2},\qquad \chi=-\frac{\Delta}{12}
$$

を得る。各 topology の特殊化値を Table 1 に示す。

**Table 1.** Topology conventions, $\chi$ values, and caveats.

| topology | $(a,b,c,u,v)$ | $\chi$ | reduction status | caveat |
|----------|--------------|-------:|-----------------|--------|
| $S^3$ | $(4,4,4,0,0)$ | $4$ | closed isotropic reference | — |
| $T^3$ | $(0,0,0,0,0)$ | $0$ | flat reference | — |
| $Nil^3$ | $(0,0,-1,0,0)$ | $-1/12$ | isotropic-scale local branch | anisotropic EOM は将来課題 |
| $Sol^3$ | $(0,0,0,1,-1)$ | $-1/3$ | isotropic-scale local branch | compact quotient は本稿範囲外 |

四つの geometry を共通の reduced-topology table のもとで扱う。 $Nil^3$ および $Sol^3$ については、本稿で用いる isotropic-scale reduction の範囲における entry を示しており、global quotient と anisotropic 拡張は今後の課題として残す。

### 2.3 EC+NY reduced action

Lorentzian EC+NY reduced action [7,8] を topology $\Sigma$ と torsion mode ごとに

$$
S_\Sigma = \int dt\, A_\Sigma\, q\!\left(N\,F_\Sigma(q,\eta,V;\chi) - \frac{\dot q^2}{N}\right)
$$

の形で書く。 $A_\Sigma$ は spatial normalization factor であり $\chi$ に依存しない係数として扱う。 $F_\Sigma$ は topology と torsion branch によって決まる effective potential であり、その具体形は Section 3 で用いる。補助場 $\eta(t)$, $V(t)$ は action 中に速度項を持たず、auxiliary sector として代数的に処理する。

Nieh-Yan term $\theta_{\rm NY}$ [8] は AX/VT では bulk $F_\Sigma$ に現れず、MX では mixing term $\propto V\eta\,\kappa^2\theta_{\rm NY}$ としてのみ現れる。本稿では bulk dynamical 解析を primary とし、 $Q_{\rm NY}$ を endpoint / boundary channel として分離する（Section 4.4）。
