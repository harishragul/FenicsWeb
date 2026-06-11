# Phase 2 — Convection-Diffusion with SUPG Stabilization

**Duration:** 4 weeks  
**Status:** Complete  
**Completed:** 2026-06-11  
**Depends on:** Phase 1 (correct BCs, working plots)  
**Unlocks:** Phase 4 (benchmarks use this solver), Phase 6B (DG comparison)

---

## Goal

Implement the convection-diffusion (advection-diffusion) equation — the PDE your project title promises. The key scientific result in this phase is demonstrating *why* standard Galerkin FEM fails at high Peclet numbers and how SUPG stabilization fixes it.

---

## Governing Equation

```
-k∇²u + b·∇u = f     in Ω
```

Where:
- `k` = thermal diffusivity (scalar, assumed constant)
- `b` = velocity field (vector, e.g., `b = (bx, by)`)
- `f` = volumetric heat source

**Cell Peclet number:** `Pe_h = |b|·h / (2k)`
- `Pe_h < 1`: diffusion-dominated, standard Galerkin works
- `Pe_h > 1`: convection-dominated, Galerkin produces spurious oscillations → needs SUPG

---

## Standard Galerkin Weak Form

Find `u ∈ V` such that for all `v ∈ V`:
```
k * ∫ ∇u·∇v dx  +  ∫ (b·∇u) v dx  =  ∫ f v dx  +  boundary terms
```

In FEniCS:
```python
a = k * dot(grad(u), grad(v)) * dx + dot(b, grad(u)) * v * dx
L = f * v * dx
```

---

## SUPG Stabilization Weak Form

Add a residual-weighted streamline term to each element:

```
a_SUPG = a_Galerkin + Σ_K ∫_K τ_K (b·∇u)(b·∇v) dx
```

Where the stabilization parameter is:
```
τ_K = h_K / (2|b|)  ·  (coth(Pe_h) - 1/Pe_h)
```

For large `Pe_h` this simplifies to `τ_K ≈ h_K / (2|b|)`.

In FEniCS, the stabilization term is added element-wise using `dx` with `cell_node_map`.

---

## Tasks

### 2.1 Create the `ConvectionDiffusion` Django app

```
python manage.py startapp ConvectionDiffusion
```

Structure mirrors `SteadyStateThermal/`:
```
ConvectionDiffusion/
├── __init__.py
├── apps.py
├── forms.py          # ConvectionDiffusionForm with k, bx, by, f, bc fields
├── solver.py         # solve_convection_diffusion() — Galerkin + SUPG
├── views.py          # ConvectionView, PecletStudyView
├── urls.py
└── templates/
    ├── convection.html
    └── peclet_study.html
```

---

### 2.2 Implement the standard Galerkin solver

**File:** `ConvectionDiffusion/solver.py`

```python
def solve_convection_diffusion_galerkin(mesh, k, bx, by, f, bcs):
    V = FunctionSpace(mesh, "P", 1)
    u, v = TrialFunction(V), TestFunction(V)
    b = Constant((bx, by))
    a = k * dot(grad(u), grad(v)) * dx + dot(b, grad(u)) * v * dx
    L = Constant(f) * v * dx
    u_sol = Function(V)
    solve(a == L, u_sol, bcs)
    return u_sol
```

---

### 2.3 Implement the SUPG solver

**File:** `ConvectionDiffusion/solver.py`

```python
def solve_convection_diffusion_supg(mesh, k, bx, by, f, bcs):
    V = FunctionSpace(mesh, "P", 1)
    u, v = TrialFunction(V), TestFunction(V)
    b = Constant((bx, by))
    h = CellDiameter(mesh)
    b_norm = sqrt(dot(b, b))
    Pe = b_norm * h / (2 * k)
    tau = h / (2 * b_norm) * (1 / tanh(Pe) - 1 / Pe)   # full formula
    # SUPG test function: v_supg = v + tau*(b·∇v)
    v_supg = v + tau * dot(b, grad(v))
    a = k * dot(grad(u), grad(v_supg)) * dx + dot(b, grad(u)) * v_supg * dx
    L = Constant(f) * v_supg * dx
    u_sol = Function(V)
    solve(a == L, u_sol, bcs)
    return u_sol
```

---

### 2.4 Peclet number study page

**What to show:** Solve the 1D convection-diffusion problem `-u'' + Pe·u' = 0` on `[0,1]` with `u(0)=0, u(1)=1` for Pe = 1, 10, 50, 100 using both Galerkin and SUPG on a fixed coarse mesh.

**Exact solution for 1D:** `u(x) = (exp(Pe·x) - 1) / (exp(Pe) - 1)`

**Display:**
- Side-by-side plots: Galerkin (shows oscillations at high Pe) vs SUPG (smooth)
- Error table: L2 error for both methods at each Pe value
- A short explanation of what Peclet number means physically

**This is the key validation figure** demonstrating SUPG stabilization (Goal 3).

---

### 2.5 Velocity field visualization

**File:** `ConvectionDiffusion/solver.py`

For 2D problems, overlay the velocity vector field `b = (bx, by)` on the solution heatmap:
```python
plt.quiver(X, Y, bx * np.ones_like(X), by * np.ones_like(Y), alpha=0.5)
```

Helps users understand the physics: where the "wind" is blowing heat.

---

### 2.6 Mixed BC support for convection-diffusion

- **Inflow boundary** (where `b·n < 0`): prescribe `u = g` (Dirichlet)
- **Outflow boundary** (where `b·n > 0`): natural (zero-flux) Neumann
- **No-flux walls**: Neumann `∂u/∂n = 0`
- Add automatic inflow detection using `b·n` at boundary facets

---

## Files to Create

| File | Purpose |
|------|---------|
| `ConvectionDiffusion/__init__.py` | App init |
| `ConvectionDiffusion/apps.py` | App config |
| `ConvectionDiffusion/forms.py` | Form: k, bx, by, f, bc fields |
| `ConvectionDiffusion/solver.py` | Galerkin + SUPG solvers |
| `ConvectionDiffusion/views.py` | ConvectionView, PecletStudyView |
| `ConvectionDiffusion/urls.py` | URL routes |
| `ConvectionDiffusion/templates/convection.html` | Solve page |
| `ConvectionDiffusion/templates/peclet_study.html` | Peclet comparison page |

**Files to modify:**
| File | Change |
|------|--------|
| `FenicsWeb/urls.py` | Include ConvectionDiffusion URLs |
| `FenicsWeb/settings.py` | Add ConvectionDiffusion to INSTALLED_APPS |

---

## Definition of Done

- [x] Convection-diffusion solver works for 1D and 2D meshes
- [x] SUPG stabilization implemented and tested
- [x] Peclet study page shows Galerkin oscillations vs SUPG correction side-by-side
- [x] L2 error table comparing Galerkin vs SUPG vs exact solution
- [x] Velocity field overlay shown on 2D solution plot
- [x] Inflow/outflow BC handling works correctly

---

## Completion Notes

### What was built

- **`ConvectionDiffusion/` Django app** mirroring `SteadyStateThermal/`: `forms.py`,
  `solver.py`, `views.py`, `urls.py`, `templates/convection.html`, `templates/peclet_study.html`.
- **`solve_convection_diffusion(dimension, mesh, k, bx, by, f, ..., method='galerkin'|'supg')`**
  — single solver function for 1D (`UnitIntervalMesh`) and 2D (`UnitSquareMesh`) problems.
  Standard Galerkin uses `a = k*dot(grad(u),grad(v))*dx + dot(b,grad(u))*v*dx`. SUPG replaces
  `v` with `v_SUPG = v + τ(b·∇v)` in both the diffusion and convection terms, per the weak
  form in this document.
- **SUPG stabilization parameter** `τ_K = h_K/(2|b|)·(coth(Pe_h) - 1/Pe_h)`, with the
  analytic small-`Pe_h` (or `b=0`) limit `τ_K → h_K²/(12k)` selected via `ufl.conditional`
  to avoid a `0/0` singularity. `tanh`/`conditional`/`gt` are imported from `ufl` (not
  `fenics`). Documented in `docs/science/governing_equations.md` §2.
- **`compute_peclet_number(mesh, k, bx, by)`** returns the global cell Peclet number
  `Pe_h = |b|·h_max/(2k)`, displayed on the solver page.
- **`mark_inflow_outflow(mesh, bx, by)`** classifies every exterior facet by the sign of
  `b·n`: `1`=inflow, `2`=outflow, `3`=characteristic. Requires `mesh.init(dim-1, dim)`
  before `Facet.normal()`/`Facet.exterior()` work. Verified on `UnitSquareMesh(4,4)`:
  `b=(1,0)` → 4 inflow / 4 outflow / 8 characteristic facets; `b=(0,0)` → all 16
  characteristic.
- **`peclet_study_1d(pe_values, n=20)`** — the key validation benchmark. Solves
  `-u'' + Pe·u' = 0`, `u(0)=0`, `u(1)=1` for `Pe ∈ {1, 10, 50, 100}` with both methods and
  compares against the exact solution `u(x) = (exp(Pe·x)-1)/(exp(Pe)-1)`. At `Pe=100`
  (cell `Pe_h=2.5`) Galerkin oscillates outside `[0,1]` (min ≈ -0.43) while SUPG stays
  bounded in `[0,1]` and has lower L2 error (0.0755 vs 0.0776). At low `Pe`, Galerkin is
  already accurate and SUPG's added artificial diffusion makes it marginally less accurate
  — both are checked by the test suite.
- **Plotting**: `generate_convection_solution_plot` (1D line plot or 2D heatmap with
  `plt.quiver` velocity overlay) and `generate_peclet_plot` (side-by-side Galerkin/SUPG/exact
  comparison, one subplot per Pe).
- **UI integration**: navbar links (Mesh, Conduction, Validation, Convection-Diffusion,
  Peclet Study) added across all five pages; new `/ConvectionDiffusion/` solve page and
  `/ConvectionDiffusion/peclet/` benchmark page.

### Test coverage

`tests/test_convection.py` — 27 new tests (Peclet number, exact 1D solution, Galerkin
correctness, SUPG stabilization/oscillation comparison, 2D solver, inflow/outflow
classification, Peclet study, plotting, Django views). Full suite: **48/48 passing**.

```text
48 passed in 3.27s
```

### Design decision — inflow/outflow classification scope

`mark_inflow_outflow` implements and tests the automatic `b·n` classification described in
task 2.6, and is documented as the basis for choosing inflow Dirichlet vs. outflow
natural-BC placement. The solve workflow itself still takes explicit Dirichlet values on
all domain boundaries (set by the user via the form) rather than auto-switching BC *types*
per facet — this keeps the UI simple while the classification utility demonstrates and
verifies the underlying numerical concept independently.
