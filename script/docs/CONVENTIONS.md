# CONVENTIONS — dppu Engine Geometric Conventions

This document fixes the **geometric, index, and sign conventions** shared between `BaseFrameEngine` and each topology runner in the `dppu` package.
**All runners must define `metric_frame` and `structure_constants` according to these conventions.**

⇒ [日本語版](CONVENTIONS_ja.md) | [SymPy ガイドライン](SymPy_guideline.md)

## 1. Scope and Assumptions

* All quantities treated here are expressed as components on a **Frame (Orthonormal Basis)**.
* `metric_frame` is the frame metric $(g_{ab})$: Euclidean paths use $(g_{ab}=\delta_{ab})$ (`Matrix.eye(dim)`), while the Lorentzian/LZ-native path uses $\eta_{ab}=\mathrm{diag}(-1,+1,+1,+1)$.
* The normal `BaseFrameEngine` pipeline primarily targets settings where **frame directional derivatives are unnecessary**.
  Specifically, runners adopt settings where structure constants $(C^{a}{}\_{bc})$ and connection coefficients $(\Gamma^{a}{}\_{bc})$ are treated as constants with respect to the frame, such as left-invariant frames.
  * **Note:** the low-level function `compute_riemann_tensor(..., frame_deriv=...)` accepts frame directional derivatives for local spin-2 jet analyses. The normal pipeline does not pass this hook.

## 2. Indices and Array Index Order

Array storage is fixed as follows:

* Structure constants: `C[a,b,c] = C^a_{bc}`
* Connection coefficients: `Gamma[a,b,c] = Γ^a_{bc}`
* Riemann curvature: `Riemann[a,b,c,d] = R^a_{bcd}`

Meaning of indices:

* $(a)$: Output (superscript) component
* $(b,c,d)$: Input (subscript) components
  In particular, $(\Gamma^{a}{}\_{bc})$ corresponds to $(\nabla\_{E_c} E_b = \Gamma^{a}{}\_{bc} E_a)$ (the last index $(c)$ is the "direction of differentiation").

## 2.5. Lorentzian/LZ-native Signature Layer

Lorentzian reductions and local jet analyses use the following standard settings:

| Item | Convention |
|---|---|
| `signature` | `"lorentzian"` |
| `frame_convention` | `"lz_native"` |
| Frame metric | $\eta_{ab}=\mathrm{diag}(-1,+1,+1,+1)$ |
| Time index | `0` |
| Spatial indices | `{1, 2, 3}` |

In this path, raising and lowering indices must not be treated as identity operations. Sign-sensitive quantities such as contortion, the Weyl scalar, and Pontryagin/Hodge contractions must use `metric` / `metric_inv` explicitly.

Implementation mapping:

* `make_config(..., signature="lorentzian")` chooses `frame_convention="lz_native"` by default.
* `dppu/topology/lz_adapter.py` maps legacy arrays into the LZ-native index order.
* Lorentzian Pontryagin uses the metric-aware implementation in `dppu/curvature/pontryagin_lz.py`, with `epsilon_tensor_up/down` rather than `epsilon_4d`.
* `hodge_dual_2form` in `dppu/curvature/hodge.py` computes the 2-form Hodge dual with metric awareness.

## 3. Definitions of Frame, Coframe, and Structure Constants (Most Important)

Let the dual of the frame be $(\{E_a\})$ and the coframe (1-forms) be $(\{e^a\})$.

### 3.1 Structure Equations of Coframe (Fixed)

$$
de^a = \frac12 C^a{}_{bc} e^b\wedge e^c,
\qquad C^a{}_{bc} = - C^a{}_{cb}.
$$

### 3.2 Commutation Relations of Dual Frame (Equivalent)

The above definition is equivalent to:

$$
[E_b, E_c] = - C^{a}{}_{bc} E_a.
$$

> Note: Many textbooks adopt $([E_b,E_c]=+f^{a}{}\_{bc}E_a)$.
> This project adopts the convention **$(C^{a}{}\_{bc}=-f^{a}{}\_{bc})$ relative to that $(f^{a}{}\_{bc})$**.

### 3.3 Runner Implementation Rules (Recommended)

* **Do not manually input C**. If possible, specify `de^a` on the runner side, extract `C^a_{bc}` by coefficient comparison, and put it into `self.data['structure_constants']=C`.
* At a minimum, automatically check that $C^a_{bc}$ is **antisymmetric in b,c**.

## 4. Connection (Spin Connection) and Metric Compatibility

Define the connection 1-form as:

$$
\omega^{a}{}\_b = \Gamma^{a}{}\_{bc} e^c
$$

Metric compatibility (Lorentz connection / Orthogonal connection) is fixed as a specification:

$$
\omega_{ab} = -\omega_{ba}
\quad(\Leftrightarrow\quad
\Gamma_{abc} = -\Gamma_{bac})
$$

where $(\Gamma_{abc} = g_{ad}\Gamma^d{}_{bc})$.

## 5. Levi-Civita Connection (General Koszul Implementation in the Current Engine)

When the frame is orthonormal and the above structure constant conventions are adopted, the Levi-Civita connection is calculated in the current engine using the following **General Koszul Formula**:

$$
\Gamma^a{}_{bc}
= \frac12\Big(
C^a{}_{bc} + C^c{}_{ba} - C^b{}_{ac}
\Big).
$$

(This is the form when the sign of the commutation relation in Section 3.2 is adopted.)

**Important Notes:**

1. This formula **does not assume a bi-invariant metric**.
   It functions correctly as the Levi-Civita connection on left-invariant frames, even for non-bi-invariant cases like Nil³.

2. In special cases where structure constants are totally antisymmetric in lower indices $C_{abc} = -C_{bac} = -C_{acb}$ (like SU(2)),
   this formula reduces to $\Gamma^a_{bc} = \frac{1}{2} C^a_{bc}$.

3. The engine automatically verifies **metric compatibility** $\Gamma_{abc} + \Gamma_{bac} = 0$ after calculation.
   If this is violated, it immediately throws an exception as an implementation error.

## 6. Torsion and Curvature

Torsion 2-form:

$$
T^a = de^a + \omega^{a}{}\_b \wedge e^b,
\qquad
T^a = \frac12 T^{a}{}\_{bc} e^b\wedge e^c.
$$

Curvature 2-form:

$$
R^{a}{}\_b = d\omega^a{}_b + \omega^a{}_c\wedge \omega^c{}_b,
\qquad
R^{a}{}\_b = \frac12 R^{a}{}\_{bcd} e^c\wedge e^d.
$$

## 7. Curvature Component Formula (Form used by Engine)

The normal pipeline uses the following form for $R^a_{bcd}$:

$$
R^{a}{}\_{bcd} = \Gamma^{a}{}\_{ec}\Gamma^{e}{}\_{bd}-\Gamma^{a}{}\_{ed}\Gamma^{e}{}\_{bc}+\Gamma^{a}{}\_{be} C^{e}{}\_{cd}.
$$

> Important: Generally, a frame directional derivative term
> $(E_c(\Gamma^{a}{}\_{bd}) - E_d(\Gamma^{a}{}\_{bc}))$
> appears here. The normal `BaseFrameEngine` pipeline does not pass this term,
> so runners should adopt settings where **$(\Gamma)$ is treated as constant in the frame direction**, such as left-invariant frames.
> local jet scripts that need time-direction derivatives call the low-level `compute_riemann_tensor(..., frame_deriv=...)` API directly.

## 8. Mandatory Self-Checks (Consistency the Runner Must Satisfy)

The runner must satisfy the following (failure implies definition inconsistency with the engine):

1. Antisymmetry of Structure Constants:

$$
C^{a}{}\_{bc} + C^{a}{}\_{cb} = 0.
$$

2. Metric Compatibility (Orthogonal Connection):

$$
\omega_{ab} + \omega_{ba} = 0.
$$

3. Antisymmetry of Riemann (Target of engine's strict check):

$$
R_{ab cd} = -R_{ba cd},\qquad
R_{ab cd} = -R_{ab dc}.
$$

---

## 9. Weyl Tensor and Conformal Scalar

### 9.1 Weyl Tensor Definition (4D)

$$
C_{abcd} = R_{abcd} - \frac{1}{2}(g_{ac}R_{bd} - g_{ad}R_{bc} - g_{bc}R_{ad} + g_{bd}R_{ac}) + \frac{R}{6}(g_{ac}g_{bd} - g_{ad}g_{bc}).
$$

Key properties:
- **Traceless**: $C^a{}_{bad} = 0$ (holds for all index pairs)
- **Conformally invariant**: invariant under $g_{ab} \to \Omega^2 g_{ab}$
- **Conformal flatness criterion**: $C_{abcd} = 0 \Leftrightarrow$ conformally flat

Note on frame basis (orthonormal): since $g_{ab} = \delta_{ab}$, raising and lowering indices is the identity operation; components $C_{abcd}$ and $C^{abcd}$ can be treated identically.

### 9.2 Weyl Scalar

$$
C^2 = C_{abcd}\,C^{abcd} = \sum_{a,b,c,d} C_{abcd}^2.
$$

In the frame basis, the sum is computed directly without index raising.
$C^2 = 0$ holds for isotropic $S^3 \times S^1$ at $\varepsilon = 0$ and for $T^3 \times S^1$ identically.

### 9.3 Engine Implementation

- Module: `dppu/curvature/weyl.py`
  - `compute_weyl_tensor(R_abcd, Ricci, R_scalar, metric, dim)` → $C_{abcd}$
  - `compute_weyl_scalar(C_abcd, metric_inv, dim)` → $C^2$
- Pipeline step: `E4.3b` (immediately after the Levi-Civita curvature computation)

---

## 10. Squashed Homogeneous Spaces and the $\varepsilon$-Parameter

### 10.1 Definition of Squashing

A volume-preserving anisotropy deformation parameter $\varepsilon$ is introduced.
$\varepsilon = 0$ is the isotropic reference point. Physical range: $\varepsilon \in (-1, +\infty)$ ($\varepsilon = -1$ is a singularity where structure constants diverge).

### 10.2 Topology-Dependent Structure Constant Scaling

The base structure constants $C^a{}_{bc}(\varepsilon=0)$ on the left-invariant frame are scaled as follows.

**$S^3 \times S^1$ (SU(2)):**

| Frame index $a$ | Scale factor $\lambda_a(\varepsilon)$ |
|---|---|
| $a \in \{0, 1\}$ | $(1+\varepsilon)^{2/3}$ |
| $a = 2$ | $(1+\varepsilon)^{-4/3}$ |

Volume preservation: $\lambda_0\lambda_1\lambda_2 = (1+\varepsilon)^{2/3+2/3-4/3} = 1$.

**$\mathrm{Nil}^3 \times S^1$ (Heisenberg group):**

| Non-trivial index $a$ | Scale factor |
|---|---|
| $a = 2$ (non-commutative component) | $(1+\varepsilon)^{-4/3}$ |

As $\varepsilon \to +\infty$, the structure constants vanish and the geometry asymptotes to flat $T^3$-like behaviour.

**$T^3 \times S^1$ (Abelian group):**

Structure constants are identically zero ($C^a{}_{bc} = 0$). The $\varepsilon$-deformation is not defined. $C^2 = 0$ holds identically, so the Weyl term vanishes (null test).

### 10.3 Physical Limits

| Limit | Physical meaning |
|---|---|
| $\varepsilon = 0$ | Isotropic $S^3$: $C^2 = 0$, stable vacuum of Paper I |
| $\varepsilon \to +\infty$ | $\mathrm{Nil}^3$ flat limit ($C^2 \to 0$ asymptotically) |
| $\varepsilon \to -1^+$ | Singularity of $S^3$ structure (physically inaccessible) |
| $\varepsilon = -2$ | Mathematical root of $C^2 = 0$, excluded since $\varepsilon < -1$ |

> $\varepsilon$ is the coordinate along $T_1$ in the full spin-2 quintet (see §11). Phase transitions occur at $\varepsilon_{c+} \approx 0.483$ (squash softening) and $\varepsilon_{c-} \approx -0.293$ (shear condensation).

---

## 11. Full Shear Tensor Basis and Spin-2 Quintet

### 11.1 Traceless Symmetric Tensor Basis $T_A$

The full shear deformation of $S^3$ is parameterised by a traceless symmetric $3\times3$ matrix $h_{ij}$ (spin-2 quintet). The orthonormal basis $\{T_A\}\_{A=1}^{5}$ satisfies $\mathrm{Tr}(T_A T_B) = \delta\_{AB}$:

| Index $A$ | Matrix | SO(2) charge (squash axis = 2-axis) |
|---|---|---|
| 1 | $\mathrm{diag}(1,1,-2)/\sqrt{6}$ | $m=0$ (singlet) |
| 2 | $\mathrm{diag}(1,-1,0)/\sqrt{2}$ | $m=\pm 2$ doublet |
| 3 | $(e_0 e_1^T + e_1 e_0^T)/\sqrt{2}$ | $m=\pm 2$ doublet |
| 4 | $(e_0 e_2^T + e_2 e_0^T)/\sqrt{2}$ | $m=\pm 1$ doublet |
| 5 | $(e_1 e_2^T + e_2 e_1^T)/\sqrt{2}$ | $m=\pm 1$ doublet |

In code (SymPy):
```python
T[1] = Matrix([[1, 0, 0], [0, 1, 0], [0, 0, -2]]) / sqrt(6)
T[2] = Matrix([[1, 0, 0], [0, -1, 0], [0, 0, 0]]) / sqrt(2)
T[3] = Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 0]]) / sqrt(2)
T[4] = Matrix([[0, 0, 1], [0, 0, 0], [1, 0, 0]]) / sqrt(2)
T[5] = Matrix([[0, 0, 0], [0, 0, 1], [0, 1, 0]]) / sqrt(2)
```

### 11.2 Embedding of $(\varepsilon, s)$ in the Quintet

The squash parameter $\varepsilon$ and shear parameter $s$ map to coordinates along $T_1$ and $T_2$:

$$
h = q_1 T_1 + q_2 T_2 + \cdots,\quad
q_1 = \tfrac{\sqrt{6}}{2}\ln(1+\varepsilon),\quad
q_2 = \sqrt{2}\ln(1+s).
$$

The two-parameter frame used by `S3Geometry` is:
```
factor_0 = (1+ε)^{2/3} (1+s)²
factor_1 = (1+ε)^{2/3} / (1+s)²
factor_2 = (1+ε)^{−4/3}             (independent of s)
```

### 11.3 Isotropic Point and Schur Degeneracy

At the isotropic point ($\varepsilon = s = 0$), SO(3) symmetry forces the Hessian to be proportional to the identity:

$$
H_{AB} = \mu_q^2\,\delta_{AB},\quad \mu_q^2 = 48\pi^2.
$$

The off-diagonal Hessians $H_{AB}$ ($A \neq B$) and the cross-Hessians between spin-2 and other sectors (spin-0, spin-1) all vanish at the isotropic point.

### 11.4 SO(3) → SO(2) Splitting

For $\varepsilon \neq 0$ the symmetry reduces to SO(2), splitting the 5-fold degeneracy into at most 3 distinct eigenvalues:

| Sector | Critical point | Physical meaning |
|---|---|---|
| $m=0$ singlet ($T_1$) | $\varepsilon_{c+} \approx 0.483$ | squash softening (1st soft mode) |
| $m=\pm2$ doublet ($T_2, T_3$) | $\varepsilon_{c-} \approx -0.293$ | 2nd-order transition to $s \neq 0$ |
| $m=\pm1$ doublet ($T_4, T_5$) | massive at both $\varepsilon_{c\pm}$ | no transition |

---

## 12. Spin-1 Sector: Physical/Auxiliary Decomposition

### 12.1 Raw Variables

The spin-1 sector has 6 raw DOFs: twist $\omega_k$ and mixing $\delta_k$ ($k=0,1,2$). These correspond to `omega1, omega2, omega3` and `delta0, delta1, delta2` in the engine (BOTH fiber mode).

### 12.2 Null Mode and Non-Dynamical Auxiliary

The field-space metric $G$ (kinetic metric, rank-3 in 6D) has a null direction:

$$
Y_k \propto \left(\omega_k,\, \delta_k\right) = \left(1,\, \frac{2L}{r_0}\right).
$$

**Classification — Case B (confirmed)**: $Y_k$ is a non-dynamical auxiliary direction (primary constraint analog). It is **not** a gauge symmetry and **not** a variable redundancy ($\delta g_{\mu\nu} \neq 0$ along $Y_k$).

### 12.3 Physical Triplet $X_k$

The unique propagating spin-1 triplet:

$$
X_k = \frac{-2\omega_k + 3\delta_k}{\sqrt{13}},
$$

with generalized eigenvalue $\lambda_{\rm phys} = H_{\rm phys}/G_{\rm phys} \approx 4.19 > 0$ (3-fold degenerate, SO(3)-symmetric).

### 12.4 Generalized Eigenvalue Problem

Stability of the spin-1 sector is assessed via:

$$
H\,v = \lambda\,G\,v.
$$

Procedure:
1. Compute $H$ (6×6 Hessian of $V_{\rm eff}$ with respect to $(\omega_k, \delta_k)$ ) by finite differences.
2. Compute $G$ (6×6 field-space metric) via the velocity-mode engine (`enable_velocity=True`).
3. Identify the null direction $Y_k$ of $G$; project $H$ and $G$ onto the orthogonal complement.
4. Solve $H_{\rm phys} v = \lambda_{\rm phys} G_{\rm phys} v$ in the 3D physical subspace.

> **Note**: The negative eigenvalue of $H$ in the full 6×6 problem lies in the null direction of $G$ and is non-propagating. It does **not** indicate a tachyon.

---

## 13. UnifiedEngine and DOFConfig Architecture

### 13.1 Class Hierarchy

The current engine consolidates topology engines into a unified interface:

```
BaseFrameEngine          (dppu/engine/pipeline.py)
    └── TopologyEngine   (dppu/topology/base_topology.py)  ← abstract
            ├── S3Geometry    (dppu/topology/s3.py)
            ├── T3Geometry    (dppu/topology/t3.py)
            ├── Nil3Geometry  (dppu/topology/nil3.py)
            └── Sol3Geometry  (dppu/topology/sol3.py)
```

`UnifiedEngine` (in `dppu/topology/unified.py`) is a factory that selects the appropriate subclass based on a `DOFConfig`.

### 13.2 DOFConfig Parameters

| Parameter | Type | Meaning |
|---|---|---|
| `topology` | `TopologyType` | `S3` / `T3` / `NIL3` / `SOL3` |
| `enable_squash` | `bool` | Activate $\varepsilon$ deformation |
| `enable_shear` | `bool` | Activate $s$ (T₂) shear |
| `enable_offdiag_shear` | `bool` | Activate $q_3,q_4,q_5$ off-diagonal shear |
| `fiber_mode` | `FiberMode` | `NONE` / `TWIST` / `MIXING` / `BOTH` |
| `isotropic_twist` | `bool` | Use single $\omega$ for all three directions |
| `torsion_mode` | `Mode` | `AX` / `VT` / `MX` |
| `ny_variant` | `NyVariant` | Nieh-Yan variant |
| `enable_velocity` | `bool` | Activate velocity symbols for $G$ metric computation |
| `skip_antisymmetry_check` | `bool` | Skip expensive Riemann antisymmetry verification |
| `weyl_source` | `CurvatureSource` | Use `LC` or `EC` Weyl scalar in the action |
| `riemann_source` | `CurvatureSource` | Default Riemann source for diagnostics/lambdify helpers |
| `signature` | `str` | `"euclidean"` / `"lorentzian"` |
| `frame_convention` | `str` | `"legacy_euclidean"` / `"lz_native"` |

### 13.3 Usage

```python
from dppu.topology import UnifiedEngine, TopologyType, FiberMode, make_config
from dppu.torsion.mode import Mode
from dppu.torsion.nieh_yan import NyVariant

cfg = make_config(
    TopologyType.S3,
    enable_squash=False,
    fiber_mode=FiberMode.BOTH,
    isotropic_twist=False,
    torsion_mode=Mode.MX,
    ny_variant=NyVariant.FULL,
    signature="lorentzian",
)
engine = UnifiedEngine(cfg)
engine.run()

Veff = engine.data['potential']
fp   = engine.get_free_params()   # dict of active SymPy Symbols
```

### 13.4 Legacy Equivalences

| Legacy class | DOFConfig equivalent |
|---|---|
| `S3S1Engine` | `S3, squash=True, shear=False, fiber=NONE` |
| `T3S1Engine` | `T3, squash=False, shear=False, fiber=NONE` |
| `Nil3S1Engine` | `NIL3, squash=True, shear=False, fiber=NONE` |

The current implementation does not provide an automatic legacy-name factory. New code should use explicit `make_config(...)` or `make_engine(...)` calls.
