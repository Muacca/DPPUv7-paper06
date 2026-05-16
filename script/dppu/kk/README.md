# KK Layer

⇒ [日本語](README_ja.md)

Module group for Kaluza-Klein photon effective theory.

## Overview

Computes the O(A²) photon effective action from 5D Einstein-Cartan gravity via two independent routes (fast shortcut and full Riemann validation). Provides coefficient extraction for Maxwell, Chern-Simons, and mass² terms across all three topologies.

**Effective action structure** (after integrating over the internal manifold):

```
S_eff ∝ ∫ d⁴x [ Σ_{i<j} k_ij F²_{ij}  +  k_CS ε^{ijk} A_i ∂_j A_k  +  Σ_i m²_i A_i² ]
```

---

## Two Computation Routes

```
Shortcut route (fast):
    field_strength → gamma_gamma → extractor

Full validation route (slow):
    field_strength → full_riemann → extractor
    └── compare with shortcut for cross-verification
```

---

## Modules

### field_strength.py

Topology-aware modified field strength F̃ and O(A¹) connection perturbation ω⁽¹⁾.

**Key functions:**

- `make_F_plain(j, k, A, dA)`: Plain field strength F_{jk} = ∂_j A_k − ∂_k A_j
- `make_F_tilde(j, k, A, dA, corrections)`: Modified field strength F̃ including topology-specific corrections
- `s3_corrections(A)`: S³ correction dict: `{(j,k,l): coeff}` (SU(2) structure constants)
- `nil3_corrections(A)`: Nil³ correction dict (Heisenberg structure constants)
- `t3_corrections(A)`: T³ correction dict (empty — flat, no corrections)
- `sol3_corrections(A, scale=1)`: Sol³ correction dict for homogeneous Sol structure constants
- `sol3_inhomogeneous_corrections(A, R_inv)`: Sol³ correction dict with an arbitrary `1/R` profile
- `make_omega1(F_tilde_fn, r0, L)`: Build O(A¹) connection perturbation via Koszul formula
- `omega1_to_array(omega1_fn, n)`: Convert ω⁽¹⁾ function to a 3D array

### gamma_gamma.py

Γ×Γ shortcut for the O(A²) Ricci scalar — the fast route.

**Key functions:**

- `gamma_gamma_ricci(omega1_fn)`: Compute R^{Γ²} = Σ_{a,b,c} ω⁽¹⁾_{abc} ω⁽¹⁾_{abc} (shortcut; exact at O(A²))

**Note:** This is algebraically equivalent to the full Riemann route at O(A²) and is ~21× faster.

### full_riemann.py

Full Riemann validation route — slow but used to cross-check the shortcut.

**Key functions:**

- `make_omega0(topology, r0, L)`: Build background (O(A⁰)) connection for the given topology
- `full_riemann_scalar_a2(omega0, omega1_fn)`: Compute O(A²) piece of full R[ω⁽⁰⁾+ω⁽¹⁾]

**Sign convention:** Returns `−R_a2` to match the `[E_b, E_c] = −C^a_{bc} E_a` convention used throughout the engine.

### extractor.py

Maxwell, Chern-Simons, and mass² coefficient extraction from the Ricci scalar.

**Key functions:**

- `extract_maxwell(R_scalar, dA)`: Extract `{(i,j): coeff}` for each F²_{ij} direction
- `extract_mass(R_scalar, A, dA)`: Extract `{i: coeff}` for each A_i² mass term
- `extract_cs(R_scalar, A, dA)`: Extract k_CS for the ε^{ijk} A_i ∂_j A_k term (or `None` if absent)
- `extract_all(R_scalar, A, dA)`: Convenience wrapper returning `{'maxwell': ..., 'mass': ..., 'cs': ...}`

**Return convention:** All extractors return dicts so anisotropic spaces (e.g. squashed S³) are handled exactly. Empty dict means the term is absent.

### validator.py

Cross-route validation utility.

**Key functions:**

- `validate_kk_routes(topology)`: Run both routes for the given topology and assert agreement

---

## Usage

### Quick Start (Shortcut Route)

```python
from sympy import symbols, Symbol
from dppu.kk import (
    s3_corrections, make_F_tilde, make_omega1,
    gamma_gamma_ricci, extract_all,
)

r0, L = symbols('r0 L', positive=True)
A   = [Symbol('A0'), Symbol('A1'), Symbol('A2')]
dA  = [[Symbol(f'dA{j}{k}') for k in range(3)] for j in range(3)]

# S³ topology
corrections = s3_corrections(A)
F_tilde_fn  = lambda j, k: make_F_tilde(j, k, A, dA, corrections)
omega1_fn   = make_omega1(F_tilde_fn, r0, L)
R_GG        = gamma_gamma_ricci(omega1_fn)
coeffs      = extract_all(R_GG, A, dA)

print(coeffs['maxwell'])   # {(0,1): ..., (0,2): ..., (1,2): ...}
print(coeffs['mass'])      # {0: ..., 1: ..., 2: ...}  (S³: 3-fold degenerate)
print(coeffs['cs'])        # k_CS expression
```

### Validation

```python
from dppu.kk.validator import validate_kk_routes

validate_kk_routes('t3')    # T³ — fast
validate_kk_routes('nil3')  # Nil³
validate_kk_routes('s3')    # S³ — slow
```

`field_strength.py` also contains Sol³ correction helpers.  The route validator currently remains scoped to `S3/T3/Nil3`.

---

## Topology Results Summary

| Topology | Maxwell k | CS k | mass² |
|----------|-----------|------|-------|
| T³×S¹ | −L²/(2r₀⁴) (uniform) | 0 | 0 (massless) |
| Nil³×S¹ | −L²/(2r₀⁴) (uniform) | −L²/r₀⁴ × A₂∂₀A₁ (01 only) | −L²/(2r₀⁴) (A₂ only) |
| S³×S¹ | −L²/(2r₀⁴) (uniform) | −2L²/r₀⁴ (all 3 dirs) | −2L²/r₀⁴ (3-fold degenerate) |

The Maxwell coefficient is the same universal value −L²/(2r₀⁴) across all three topologies.

---

## Dependencies

- [engine](../engine/README.md): Koszul formula for ω⁽¹⁾
- [topology](../topology/README.md): Background connection ω⁽⁰⁾
- [utils](../utils/README.md): Epsilon symbols (`epsilon_3d`)

## Related Usage

- Use `dppu.kk.validator.validate_kk_routes(...)` to compare the shortcut and full-Riemann routes.
- Use `dppu.kk.extractor.extract_all(...)` when a script needs Maxwell, Chern-Simons, and mass coefficients from a symbolic Ricci scalar.
