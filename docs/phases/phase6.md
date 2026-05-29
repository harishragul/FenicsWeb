# Phase 6 — Research Contribution

**Duration:** 4 weeks
**Status:** Not Started
**Depends on:** Phase 1–4 (validated solvers), Phase 5 (async infrastructure)
**Unlocks:** Phase 7 (publication — needs a novel result to report)

---

## Goal

Produce one novel scientific result that goes beyond implementing known algorithms. This is what separates a thesis from a project report and what a PhD application needs to demonstrate research potential.

The research gap identified from the literature search: **no published tool provides a web-accessible interface for inverse conductivity recovery in a convection-dominated heat transfer domain.** dolfin-adjoint exists as a library, but combining it with SUPG-stabilised convection-diffusion and sparse boundary measurements — inside a validated, documented web platform — is an original contribution.

**Choose ONE option.** Recommendation: **Option A** (sharpened with a novel angle based on the literature gap).

---

## Option A — Inverse Heat Conduction in a Convection-Dominated Domain (Recommended)

### Research Question

> *"Can we recover the unknown thermal conductivity k(x) from sparse, noisy boundary temperature measurements when the domain is convection-dominated (high Peclet number)? How does convection strength affect identifiability and regularization requirements?"*

### Why This Is Novel

Prior work on inverse HCP via adjoint methods (dolfin-adjoint) treats pure diffusion problems. The **convection-dominated setting introduces new challenges** that have not been systematically studied in the web-platform context:

- The convection term creates asymmetric information flow — measurements downstream carry less information about upstream conductivity
- SUPG stabilization of the forward problem interacts with the adjoint in a non-trivial way
- The regularization parameter `α` required for stable recovery changes with Peclet number

This identifiability-vs-Peclet study is the novel result. It does not require inventing a new algorithm — the novelty is the systematic study and the accessible web interface.

### Mathematical Formulation

**Forward problem (convection-diffusion):**

```text
-∇·(k∇u) + b·∇u = f    in Ω
```

**Inverse problem:** Given noisy measurements `u_obs = {u(xᵢ) + εᵢ}` at sparse points, find `k(x)` minimising:

```text
J(k) = Σᵢ (u(xᵢ) - u_obs_i)²  +  α ∫ |∇k|² dx
```

**Key addition over standard IHCP:** the velocity field `b` is now present. The adjoint equation becomes:

```text
-∇·(k∇p) - b·∇p = -2 Σᵢ (u(xᵢ) - u_obs_i) δ(x - xᵢ)
```

Note the sign flip on the convection term in the adjoint — this is the mathematical consequence of convection and is what makes the problem harder.

**Gradient:** `dJ/dk = ∇u · ∇p + 2α(-∇²k)` — computed automatically by dolfin-adjoint.

### Implementation Plan

**New app:** `Research/` with sub-module `Research/inverse/`

```python
# Research/inverse/solver.py
from fenics import *
from fenics_adjoint import *
from scipy.optimize import minimize

def solve_inverse_convection(mesh, b_vec, u_obs_coords, u_obs_vals,
                              f, alpha, k_init=1.0):
    """Recover k(x) from sparse measurements under convection field b."""
    V = FunctionSpace(mesh, "P", 1)
    K = FunctionSpace(mesh, "P", 1)

    k = Function(K)
    k.interpolate(Constant(k_init))
    b = Constant(b_vec)

    def forward(k):
        u = Function(V)
        v = TestFunction(V)
        # SUPG-stabilised convection-diffusion forward problem
        h = CellDiameter(mesh)
        b_norm = sqrt(dot(b, b) + 1e-12)
        Pe = b_norm * h / (2 * k)
        tau = h / (2 * b_norm) * (1 / tanh(Pe) - 1 / Pe)
        v_s = v + tau * dot(b, grad(v))
        F = (k * dot(grad(u), grad(v_s))
             + dot(b, grad(u)) * v_s
             - Constant(f) * v_s) * dx
        bc = DirichletBC(V, Constant(0), "on_boundary")
        solve(F == 0, u, bc)
        return u

    def objective(k_vals):
        k.vector()[:] = k_vals
        u = forward(k)
        misfit = sum((u(xi) - ui) ** 2
                     for xi, ui in zip(u_obs_coords, u_obs_vals))
        reg = alpha * assemble(dot(grad(k), grad(k)) * dx)
        return float(misfit + reg)

    result = minimize(objective, k.vector()[:], method='L-BFGS-B',
                      options={'maxiter': 200})
    k.vector()[:] = result.x
    return k
```

### Study Design (the novel part)

#### Experiment 1 — Identifiability vs Peclet number

Fix `k_true(x) = 1 + 0.5 sin(2πx)`, 10 measurement points, σ = 0.01 noise.

Vary Pe = 0 (pure diffusion), 1, 5, 10, 50 by changing `|b|`.

For each Pe, find the optimal `α` using the L-curve and report:

- L2 error in k recovery
- Number of optimizer iterations to convergence
- Whether the recovery degrades upstream or downstream of the flow

Expected result: recovery quality degrades with Pe; downstream measurements become uninformative; a larger `α` is needed to stabilize recovery at high Pe.

#### Experiment 2 — Measurement sparsity study

Fix Pe = 10. Vary number of measurements: 3, 5, 10, 20 points.

Report L2(k) error vs number of measurements. Find the minimum number of measurements required for stable recovery at this Pe.

#### Experiment 3 — Noise robustness

Fix Pe = 10, 10 measurements. Vary noise σ = 0, 0.001, 0.01, 0.05, 0.1.

Report L2(k) error vs σ at the optimal `α` for each noise level.

### Results to Report

- Figure: `k_true` vs `k_recovered` for Pe = 0, 5, 50 side-by-side
- Figure: L-curve for Pe = 10 (data misfit vs regularization norm)
- Table: L2(k) error vs Pe at optimal α
- Table: L2(k) error vs number of measurement points
- Table: L2(k) error vs noise level σ
- Discussion: physical interpretation — why does high Pe hurt recovery?

These results are original. They can form a standalone section in the thesis and strengthen the JOSS paper's claim of research use.

---

## Option B — SUPG vs Discontinuous Galerkin Comparison

### Option B Research Question

> *"For convection-dominated heat transfer, how do SUPG-CG and Interior Penalty DG methods compare in accuracy, stability, and computational cost across Peclet numbers?"*

### Option B Publishability

A direct, quantitative comparison with convergence tables and CPU timing is a contribution in itself. The JOSS paper can cite this study as evidence of research use.

### DG Weak Form (Interior Penalty)

```python
V = FunctionSpace(mesh, "DG", 1)
n = FacetNormal(mesh)
h_avg = (CellDiameter(mesh)('+') + CellDiameter(mesh)('-')) / 2
sigma = Constant(10.0)   # penalty parameter
un = 0.5 * (dot(b, n) + abs(dot(b, n)))

a = (k * dot(grad(u), grad(v)) * dx
     - k * dot(n, avg(grad(u))) * jump(v) * dS
     - k * dot(n, avg(grad(v))) * jump(u) * dS
     + (sigma / h_avg) * jump(u) * jump(v) * dS
     + dot(b, grad(u)) * v * dx
     - (un('+') - un('-')) * jump(v) * dS)
```

### Comparison Metrics

- L2 and H1 convergence rates under h-refinement
- Solution profiles at Pe = 1, 10, 100 on coarse mesh (qualitative stability)
- CPU time per solve at equivalent accuracy level
- Degrees of freedom for same L2 error

---

## Option C — Uncertainty Quantification

### Option C Research Question

> *"If thermal conductivity k is uncertain, how does that uncertainty propagate through the temperature field — and does convection amplify or dampen uncertainty?"*

### Option C Methods

Monte Carlo: Sample `k ~ N(μ, σ²)` N=1000 times, solve forward, compute `E[u]` and `Var[u]`.

Polynomial Chaos Expansion (PCE): Represent `k(ξ) = μ + σξ`, expand `u` in Hermite polynomials, solve once for all coefficients. Library: `chaospy`.

Novel angle: study how convection strength (Pe) affects the spatial distribution of `Var[u]` — does high Pe concentrate uncertainty near outflow boundaries?

---

## Phase 6 Files to Create (Option A)

| File | Purpose |
| ---- | ------- |
| `Research/__init__.py` | App init |
| `Research/apps.py` | App config |
| `Research/inverse/solver.py` | SUPG forward + adjoint inverse solver |
| `Research/inverse/forms.py` | Inverse problem parameters form |
| `Research/inverse/views.py` | InverseView, PecletStudyView, LCurveView, SparsityView |
| `Research/inverse/urls.py` | URL routes |
| `Research/inverse/templates/inverse.html` | Inverse solver page |
| `Research/inverse/templates/peclet_study.html` | Pe vs recovery quality results |
| `Research/inverse/templates/lcurve.html` | L-curve and regularization page |
| `tests/test_inverse.py` | Tests: noiseless recovery error < 1e-2 |

---

## Definition of Done (Option A)

- [ ] Forward solver (convection-diffusion with SUPG) wraps correctly with dolfin-adjoint
- [ ] Inverse solver recovers `k_true` from noiseless data with L2 error < 0.01
- [ ] Inverse solver recovers reasonable `k` from σ = 0.01 noise
- [ ] Experiment 1 complete: L2(k) error table vs Pe (5 Peclet values)
- [ ] Experiment 2 complete: L2(k) error vs number of measurement points
- [ ] Experiment 3 complete: L2(k) error vs noise level σ
- [ ] L-curve page built for Pe = 10 case
- [ ] All experiments accessible from the web UI
- [ ] `tests/test_inverse.py` passes in CI
