# SymPy Implementation Guidelines — dppu Engine

- **Target:** AI Coding Assistants & Contributors
- **Context:** Einstein-Cartan Theory on Curved Spacetime (S³×S¹, T³×S¹, Nil³×S¹)

This document establishes the implementation guidelines for symbolic computation (SymPy) and theoretical physics conventions for the `dppu` package. The following rules must be strictly adhered to.

⇒ [日本語版](SymPy_guideline_ja.md) | [Geometric Conventions](CONVENTIONS.md)

-----

## 1. Optimization Rules (Engineering)

Iron rules to prevent "computation time explosion" (taking over 2 hours) in SymPy integration calculations, ensuring completion within seconds.

### Rule 1.1: Strict use of `expand()` + `cancel()` Strategy

Do not use high-cost functions like `simplify()` on the integrand immediately before performing integration (`integrate`). Instead, perform "expansion and cancellation".

  * **Don't:**
    ```python
    density = simplify(density)  # NG: PROHIBITED. Factorizing huge expressions has extreme computational cost.
    result = integrate(density, x)
    ```
  * **Do:**
    ```python
    density = cancel(expand(density))  # OK: RECOMMENDED. Converting to a sum of polynomials induces term-by-term integration.
    result = integrate(density, x)
    ```

### Rule 1.2: Suppress Intermediate Simplification

When performing multiple integrations (e.g., $\phi$ integration $\to$ $\theta$ integration), do not apply excessive `simplify` to intermediate results. Limiting it to `cancel` and re-applying `expand` just before the final integration is faster.

### Rule 1.3: cos/sin Polynomial Auxiliary Symbols for Mixing DOF

When the mixing rotation angle $\theta_i = \arctan\delta_i$ enters the structure constants (fiber mode MIXING or BOTH), do **not** write `1/sqrt(1+delta_i**2)` directly in the pipeline. Instead, use the pre-defined **auxiliary symbols** `cd_i` and `sd_i`:

$$
cd_i \equiv \cos\theta_i = \frac{1}{\sqrt{1+\delta_i^2}},\quad
sd_i \equiv \sin\theta_i = \frac{\delta_i}{\sqrt{1+\delta_i^2}}.
$$

The engine introduces these automatically in `_build_fiber_params` for MIXING/BOTH modes. All intermediate expressions then remain **polynomials** in `{r, L, eta, V, omega_k, cd_i, sd_i, ...}`, enabling fast `expand`/`cancel` throughout the E4 pipeline.

**Why**: Repeated `1/sqrt(1+δ²)` prevents SymPy from recognizing polynomial structure, causing `cancel` to produce algebraically correct but computationally intractable forms. Replacing with `cd_i` reduces computation time from hours to minutes.

**Numerical substitution** (at $\delta=0$, i.e., isotropic mixing point):
```python
iso_subs = {fp['cd0']: 1, fp['sd0']: 0,
            fp['cd1']: 1, fp['sd1']: 0,
            fp['cd2']: 1, fp['sd2']: 0}
Veff_numeric = Veff_sym.subs(iso_subs).subs(phys_subs)
```

**Derivative with respect to $\delta_i$** at $\delta=0$: since $\partial/\partial\delta_i\big|\_0 = \partial/\partial(sd_i)\big|\_{sd=0,\,cd=1}$, differentiate by `sd_i` and then substitute `sd_i=0, cd_i=1`:
```python
dV_ddelta0 = diff(Veff_sym, fp['sd0']).subs(iso_subs)
```

### Rule 1.4: `lambdify` with CSE for Numerical Evaluation

For numerical Hessian computation and parameter scans, convert the symbolic $V_{\rm eff}$ to a fast NumPy function via `lambdify` with CSE (common subexpression elimination):

```python
from sympy import lambdify

fp = engine.get_free_params()          # ordered dict of active Symbols
sym_args = tuple(fp[k] for k in ordered_keys)

Vfunc = lambdify(sym_args, Veff_expr, modules='numpy', cse=True)
# cse=True: significant speedup for expressions with thousands of terms
```

Typical `ordered_keys` for S³ BOTH mode:
```python
ordered_keys = ['r', 'L', 'kappa', 'eta', 'V', 'theta_NY', 'alpha',
                'omega1', 'omega2', 'omega3',
                'cd0', 'sd0', 'cd1', 'sd1', 'cd2', 'sd2']
```

At the isotropic evaluation point, pass `cd_i=1, sd_i=0` for all $i$.

### Rule 1.5: Avoid `Matrix.inv()` for Large Symbolic Matrices

Calling `.inv()` on a SymPy `Matrix` containing multiple symbols can cause indefinite hangs (observed: >30 minutes for $4\times4$ matrices with polynomial entries).

* **Don't**: `G_inv = G_matrix.inv()`
* **Do**: Use an explicit closed-form inverse for small matrices, or restrict to numerical inversion after `lambdify`.

```python
# Example: explicit 2×2 inverse
det = a * d - b * c
G_inv = Matrix([[d, -b], [-c, a]]) / det
```

This issue was encountered in the mixing-frame structure constant transformation; the fix uses an explicit `G_inv` formula in `dppu/topology/s3.py`.

-----

## 2. Theoretical Implementation Rules

Iron rules for maintaining consistency in tensor operations within Curved Spacetime.

### Rule 2.1: Prohibition of Direct Index Manipulation (Robust Method)

In environments with non-diagonal metrics or where $g_{\mu\nu} \neq 1$, swapping array indices directly (e.g., `T[mu, nu, lam]`) is not equivalent to physical tensor index manipulation (raising/lowering).

  * **Don't:**
    ```python
    # NG: PROHIBITED. Risk of referencing physically incorrect components.
    term = T_tensor[mu, nu, lam]
    ```

  * **Do:** Always manipulate indices via the metric tensor $g_{\mu\nu}$.
    1.  Lower all indices to create the fully covariant form $T_{\lambda\mu\nu}$.
    2.  Perform index permutation.
    3.  Raise indices using the metric as necessary.

#### Optimization for Orthonormal Frame Basis

In an **Orthonormal Frame Basis** with identity metric, the computational overhead of raising and lowering indices may be omitted. Direct component calculation can be used for speed in the legacy Euclidean path.

However, the sign pattern $(+1, +1, -1)$ defined by the physical definition (Hehl 1976) must be strictly observed.

Do not apply this shortcut to the Lorentzian/LZ-native path. There the metric is $\mathrm{diag}(-1,+1,+1,+1)$, and contortion, the Weyl scalar, and Pontryagin/Hodge contractions must use `metric` / `metric_inv` for index operations. The snippet below is explanatory for the legacy Euclidean or identity-metric path.

```python
# ============================================================
# Golden Logic for Contortion (Frame Basis / legacy Euclidean path)
# ============================================================
# Assumption: Metric is diagonal/identity (Orthonormal Frame)
# Therefore T^a_bc and T_abc behave identically in code logic.

K_tensor = MutableDenseNDimArray.zeros(dim, dim, dim)

for a in range(dim):
    for b in range(dim):
        for c in range(dim):
            # Formula: K_abc = (1/2)(T_abc + T_bca - T_cab)
            # Note: Using T[a,b,c] directly as T_abc
            
            term = (T_tensor[a, b, c] + T_tensor[b, c, a] - T_tensor[c, a, b])
            
            val = term * Rational(1, 2)
            
            if val != 0:
                K_tensor[a, b, c] = cancel(expand(val))
```

### Rule 2.2: Consistency Check

After constructing the EC connection ($\Gamma_{\text{EC}}$), the following verification code must be executed to confirm that the mismatch count is **0**.

```python
# Torsion Consistency Check
T_verify = Gamma_EC[lam, mu, nu] - Gamma_EC[lam, nu, mu] # Hehl definition
mismatch = count(simplify(T_verify - T_original) != 0)
assert mismatch == 0
```

-----

## 3. Standard Conventions (Theoretical)

To prevent confusion during paper writing, adhere to the **Hehl (1976) Standard**.

### 3.1 Torsion Definition

Definition of the Torsion Tensor $T^\lambda_{\ \mu\nu}$:
$$T^\lambda_{\ \mu\nu} \equiv \Gamma^\lambda_{\ \mu\nu} - \Gamma^\lambda_{\ \nu\mu}$$
(Note: In the frame basis of this engine, torsion components are extracted from the torsion 2-form $T^a = de^a + \omega^{a}{}\_b\wedge e^b$ defined in Section 6 of CONVENTIONS, via the coefficient comparison $T^a = \frac{1}{2}T^{a}{}\_{bc}\,e^b\wedge e^c$.)

### 3.2 Contortion Formula

The Verified Formula for Contortion $K^\lambda_{\ \mu\nu}$ consistent with the above Torsion definition:

$$K_{\lambda\mu\nu} = \frac{1}{2} \left( T_{\lambda\mu\nu} + T_{\mu\nu\lambda} - T_{\nu\lambda\mu} \right)$$

  * Sign Pattern: **$(+1, +1, -1)$**
  * Note: Apply this formula to $T_{\lambda\mu\nu}$ (where all indices are lowered).

### 3.3 Einstein-Cartan Connection

$${\Gamma_{\text{EC}}}^\lambda_{\ \mu\nu} = {\Gamma_{\text{LC}}}^\lambda_{\ \mu\nu} + K^\lambda_{\ \mu\nu}$$

  * $\Gamma_{\text{LC}}$: Levi-Civita Connection (Christoffel symbols)
  * $K$: Contortion

-----

## 4. Torsion Ansatz and Mode Decomposition Rules

The torsion tensor on the $M^3 \times S^1$ minisuperspace ansatz is specified via three modes.

### 4.1 Mode Definitions

| Mode | Physical component | Parameters |
|---|---|---|
| `Mode.AX` | Axial component (T1) only | $\eta \neq 0$, $V = 0$ |
| `Mode.VT` | Vector-trace component (T2) only | $\eta = 0$, $V \neq 0$ |
| `Mode.MX` | Both T1 and T2 | $\eta \neq 0$, $V \neq 0$ |

### 4.2 Physical Correspondence

- **T1 (Axial component):** Dual to axial vector $S^\mu = (\eta/r)(0,0,0,1)$. For spatial indices $a,b,c \in \{0,1,2\}$: $T_{abc} = (2\eta/r)\,\varepsilon_{abc}$.
- **T2 (Vector-trace component):** Dual to vector $V_\mu = V\,\delta^3_\mu$ ($\tau$-component only). $T_{abc} = \frac{1}{3}(\delta_{ac}V_b - \delta_{ab}V_c)$.

### 4.3 Implementation Rule

Use `construct_torsion_tensor(mode, r, eta, V, metric, dim)` from `dppu/torsion/ansatz.py`.
Do **not** construct $T_{abc}$ by hand.

-----

## 5. Nieh-Yan Topological Term Variants

### 5.1 Nieh-Yan Decomposition

The full Nieh-Yan density:
$N = N_{\mathrm{TT}} - N_{\mathrm{Ree}},$
$N_{\mathrm{TT}} = \frac{1}{4}\varepsilon^{abcd}T^{e}{}\_{ab}T_{ecd},\qquad N_{\mathrm{Ree}} = \frac{1}{4}\varepsilon^{abcd}R_{abcd}.$

### 5.2 Variant Selection

| `NyVariant` | Density used |
|---|---|
| `NyVariant.TT` | $N_{\mathrm{TT}}$ only |
| `NyVariant.REE` | $N_{\mathrm{Ree}}$ only |
| `NyVariant.FULL` | $N_{\mathrm{TT}} - N_{\mathrm{Ree}}$ (canonical) |

### 5.3 Implementation

All three variants are computed in pipeline step `E4.10`. Select the variant via the `ny_variant` argument in the engine's `__init__`. See `dppu/torsion/nieh_yan.py`.

-----

## 6. Extended Lagrangian and Weyl Coupling Constant $\alpha$

### 6.1 Action Form

$$S = \int L \times \mathrm{Vol},\qquad
L = \frac{R}{2\kappa^2} + \theta_{\mathrm{NY}}\times N + \alpha\times C^2.$$

| Parameter | Meaning |
|---|---|
| $\kappa$ | Einstein-Cartan gravitational coupling |
| $\theta_{\mathrm{NY}}$ | Nieh-Yan coupling (topological) |
| $\alpha$ | Weyl coupling (conformal invariant term) |

For $\alpha \leq 0$, Theorem 1 guarantees protection of the stable vacuum. For $\alpha > 0$, Theorem 2 guarantees $V_{\rm eff} \to -\infty$.

### 6.2 Obtaining the Effective Potential

**UnifiedEngine** :
```python
engine.run()
fp   = engine.get_free_params()          # dict of active SymPy Symbols
Veff = engine.data['potential']          # raw SymPy expression

# Fast numerical evaluation (see Rule 1.4):
from sympy import lambdify
Vfunc = lambdify(tuple(fp[k] for k in ordered_keys), Veff, modules='numpy', cse=True)
```

$V_{\rm eff} = -S$. Extracted in pipeline step `E4.13`.

### 6.3 Implementation

See `compute_lagrangian()` in `dppu/action/lagrangian.py` and `dppu/action/potential.py`.

-----

## 7. Numerical Optimisation Strategy (Phase Atlas Search)

### 7.1 Two-Stage Strategy

**Stage 1: Brute-force grid search**

`scipy.optimize.brute` with `Ns` points per axis over the $(r, \varepsilon)$ 2D grid, to locate the basin of the global minimum.

**Stage 2: Multi-start L-BFGS-B refinement**

Top- $N$ candidates from Stage 1 are used as starting points for `scipy.optimize.minimize` (L-BFGS-B, `ftol=1e-8`) to obtain a high-precision minimum.

### 7.2 Stability Classification

After locating the minimum at $(r^\*, \varepsilon^\*)$ , classify as follows:

| Class | Condition | Physical meaning |
|---|---|---|
| Type-I | $V(r)$ rises as $r \to 0$ (barrier present) | Stable vacuum (nucleation barrier exists) |
| Type-II | $V(r)$ decreases monotonically as $r \to 0$ (no barrier) | Spontaneous nucleation possible |
| Type-III | No local minimum in the physical region | Unstable configuration |

Use `analyze_stability()` from `dppu/scanning/stability.py`.

### 7.3 Important Notes

- For $\alpha > 0$, the optimiser attaches to the search boundary $(r \to 0^+,\,\varepsilon \to -1^+)$ (`converged = False`). This is not an optimisation failure; it correctly reports that the potential is unbounded below within the search range.

-----

## 8. Self-Duality (SD) Diagnostic Rules

### 8.1 Hodge Dual of Curvature

$$(*R)^{ab}{}_{cd} = \frac{1}{2}\varepsilon_{cdef}\,R^{ab,ef},$$

where $\varepsilon_{cdef}$ is the Levi-Civita symbol in the frame basis (use `epsilon_4d()` from `dppu/utils/epsilon.py`).

### 8.2 Pontryagin Inner Product and SD Residuals

$$E_{RR} = \langle R, R\rangle = R_{abcd}R^{abcd},\qquad
P = \langle R, *R\rangle = R_{abcd}(*R)^{abcd}.$$

| Condition | Physical state |
|---|---|
| SD residual $< \varepsilon_{\rm SD}$ and $\|R\| > \varepsilon_R$ | Self-dual |
| ASD residual $< \varepsilon_{\rm SD}$ and $\|R\| > \varepsilon_R$ | Anti-self-dual |
| `P_form = 0` | The form-Hodge Pontryagin density is exactly zero for AX/VT/MX EC branches |
| `P_int = 0` — AX/VT | The internal-pair Pontryagin diagnostic vanishes for AX/VT |
| `P_int \neq 0` — MX | $P_{\rm int}=2V\eta(-V^2r^2+9\eta^2-36)/(9r^3)$; source is torsion mixing $V\eta$, not twist |

### 8.3 Usage

```python
from dppu.curvature.sd_extension import SDExtensionMixin
SDExtensionMixin.attach_to(engine)          # Dynamically attach methods
R = engine.get_R_ab_cd_numerical(params_dict)
diag = engine.evaluate_sd_status(params_dict)
# Check diag['P_RstarR'] for the numerical Pontryagin diagnostic value
```

See `dppu/curvature/sd_extension.py`, `dppu/curvature/sd_diagnostics.py`, and `dppu/curvature/pontryagin.py`.

-----

## 9. Numerical Hessian Analysis

### 9.1 Pattern

Assesses stability by computing numerical Hessians at the isotropic point via finite differences, rather than symbolic differentiation.

```python
import numpy as np
from sympy import lambdify

# 1. Run engine and build fast numerical function (Rule 1.4)
fp = engine.get_free_params()
Veff_sym = engine.data['potential']
Vfunc = lambdify(tuple(fp[k] for k in ordered_keys), Veff_sym, modules='numpy', cse=True)

# 2. Wrapper that evaluates V_eff at a perturbed point
phys_base = dict(r=r0, L=1, kappa=1, eta=1, V=1, theta_NY=0, alpha=0,
                 omega1=0, omega2=0, omega3=0,
                 cd0=1, sd0=0, cd1=1, sd1=0, cd2=1, sd2=0)

def veff(**overrides):
    vals = {**phys_base, **overrides}
    return float(Vfunc(*[vals[k] for k in ordered_keys]))

# 3. Central-difference Hessian
h = 1e-4
keys = ['omega1', 'omega2', 'omega3', 'sd0', 'sd1', 'sd2']  # sd_i for delta derivatives
H = np.zeros((6, 6))
for i, ki in enumerate(keys):
    for j, kj in enumerate(keys):
        H[i, j] = (veff(**{ki: +h, kj: +h}) - veff(**{ki: +h, kj: -h})
                 - veff(**{ki: -h, kj: +h}) + veff(**{ki: -h, kj: -h})) / (4 * h**2)
```

### 9.2 Derivative Convention for $\delta_i$

Because `delta_i` appears only through `cd_i` and `sd_i`, the correct finite-difference variable is `sd_i` (not `delta_i`):

$$
\frac{\partial V}{\partial \delta_i}\bigg|_{\delta=0}
= \frac{\partial V}{\partial sd_i}\bigg|_{sd=0,\,cd=1}.
$$

Pass `sd_i = ±h` (with `cd_i = 1`) in the finite-difference calls.

### 9.3 Generalized Eigenvalue Problem for Spin-1

To assess spin-1 stability, solve $H v = \lambda G v$ where $G$ is the field-space metric:

```python
from scipy.linalg import eigh

# G is computed via enable_velocity=True (see CONVENTIONS §12.4)
# G is rank-3 in 6D → project onto physical subspace first

# null direction Y_k ∝ (1, 2L/r0)  for each k in {0,1,2}
norm_Y = np.sqrt(1 + (2*L/r0)**2)
Y_blocks = [np.array([1, 0, 0, 2*L/r0, 0, 0]) / norm_Y,   # k=0
            np.array([0, 1, 0, 0, 2*L/r0, 0]) / norm_Y,   # k=1
            np.array([0, 0, 1, 0, 0, 2*L/r0]) / norm_Y]   # k=2

# Build 3×3 physical-subspace matrices via basis X_k orthogonal to Y_k

eigenvalues, _ = eigh(H_phys, G_phys)
# lambda_phys > 0 confirms no tachyon in the physical spin-1 sector
```

> **Warning**: Do not use `eigh(H_full, G_full)` directly; $G$ is singular and `scipy` will raise a `LinAlgError`. Always project to the physical subspace first.
