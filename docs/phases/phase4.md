# Phase 4 — Scientific Rigor & Benchmarks

**Duration:** 3 weeks  
**Status:** Complete  
**Depends on:** Phase 1 (correct solver), Phase 2 (convection-diffusion)  
**Unlocks:** Phase 6 (documentation cites these results), Phase 7 (publication evidence)

---

## Goal

Add the Verification and Validation (V&V) infrastructure that separates a research tool from a student project. The plots and tables produced here form the verification and validation evidence for the project (Goals 5–7).

---

## Why V&V Matters

A PhD guide or application reviewer will immediately ask: *"How do you know your solver is correct?"*

The answer requires:
1. **Verification** — does the code solve the equations correctly? (Compare to analytical solutions)
2. **Validation** — does the model represent reality? (Compare to established benchmark results from literature)

---

## Tasks

### 4.1 Create the `Validation` Django app

```
python manage.py startapp Validation
```

```
Validation/
├── __init__.py
├── apps.py
├── analytical.py     # Exact solution functions for all test problems
├── convergence.py    # h-refinement and p-refinement study runners
├── benchmarks.py     # De Vahl Davis and other benchmark problems
├── views.py          # ConvergenceView, BenchmarkView, ValidationDashboard
├── urls.py
└── templates/
    ├── convergence.html
    ├── benchmark.html
    └── dashboard.html
```

---

### 4.2 h-Refinement Convergence Study

**Purpose:** Verify that the FEM error decreases at the theoretically predicted rate as mesh is refined.

**Test problem (1D):**
- `-k u'' = f` on `[0, 1]` with `u(0)=0, u(1)=0`, `f = 2k`
- Exact: `u(x) = x(1-x)`

**Test problem (2D):**
- `-k ∇²u = f` on `[0,1]²`, Dirichlet `u=0` on all sides
- Choose `f` such that `u = sin(πx)sin(πy)` is exact
- Then `f = 2π²k · sin(πx)sin(πy)`

**Procedure:**
```python
mesh_sizes = [4, 8, 16, 32, 64]
for n in mesh_sizes:
    mesh = UnitSquareMesh(n, n)
    u_h = solve(mesh, ...)
    u_exact = interpolate(Expression("sin(pi*x[0])*sin(pi*x[1])", degree=5), V)
    L2_error = sqrt(assemble((u_h - u_exact)**2 * dx))
    H1_error = sqrt(assemble(dot(grad(u_h - u_exact), grad(u_h - u_exact)) * dx))
```

**Expected rates:**
- L2 error: `O(h²)` → slope ≈ 2.0 on log-log plot
- H1 error (energy norm): `O(h¹)` → slope ≈ 1.0 on log-log plot

**Output:** Table + log-log convergence plot (a standard verification figure in finite element analysis).

---

### 4.3 p-Refinement Study (Higher-Order Elements)

**Purpose:** Show that using higher-degree polynomials also improves accuracy, and verify super-convergence rates.

**Implementation:**
- Solve the same 2D test problem with P1, P2, P3 elements on a *fixed* mesh (n=8)
- Then run h-refinement for each element order

**Expected L2 convergence rates:**
- P1: `O(h²)`
- P2: `O(h³)`
- P3: `O(h⁴)`

**Output:** Combined convergence plot with three curves labeled P1/P2/P3 — a classic figure in numerical analysis.

**Note:** To enable P2/P3, change only one line in the solver:
```python
V = FunctionSpace(mesh, "P", element_order)   # order = 1, 2, or 3
```
Add `element_order` as a user input.

---

### 4.4 Convection-Diffusion Convergence (with SUPG)

Repeat the h-refinement study for the convection-diffusion problem (Phase 2):

**Test problem (1D):**
- `-u'' + Pe·u' = 0` on `[0,1]`, `u(0)=0, u(1)=1`
- Exact: `u(x) = (exp(Pe·x) - 1)/(exp(Pe) - 1)`

**Study:**
- Run for Pe = 0.1 (diffusion-dominated) and Pe = 100 (convection-dominated)
- Compare L2 error of Galerkin vs SUPG at each mesh size
- SUPG should converge cleanly even at high Pe; Galerkin will not on coarse meshes

---

### 4.5 De Vahl Davis Natural Convection Benchmark

**Reference:** De Vahl Davis (1983) — the standard 2D natural convection benchmark in a square cavity.

**Problem setup:**
- Square domain `[0,1]²`
- Hot left wall: `u = 1` (Dirichlet)
- Cold right wall: `u = 0` (Dirichlet)
- Insulated top and bottom: `∂u/∂n = 0` (Neumann)
- Solve the steady-state conduction problem first (the full Boussinesq flow is beyond FEniCS scope at this level, but the thermal diffusion part is directly comparable)

**Comparison metric:** Average Nusselt number on the hot wall:
```
Nu = -∫_{left wall} ∂u/∂x dy   (dimensionless heat transfer rate)
```

**Published reference values (Table 1, De Vahl Davis 1983):**
| Ra | Nu_avg |
|----|--------|
| 10³ | 1.118 |
| 10⁴ | 2.243 |
| 10⁵ | 4.519 |
| 10⁶ | 8.800 |

Compute the Nusselt number from your solver and compare in a table. Note that without the flow solver, we can only compare the thermal field at Ra=0 (pure conduction), but the infrastructure and comparison methodology is what matters for validation.

---

### 4.6 Validation Dashboard

**Template:** `Validation/templates/dashboard.html`

A single page linking to all validation studies:
- h-refinement convergence plot (1D and 2D)
- p-refinement comparison plot
- Peclet study validation table
- De Vahl Davis benchmark comparison
- Each section: result, expected value, pass/fail indicator

This dashboard is the "trust anchor" of the entire project.

---

## Files to Create

| File | Purpose |
|------|---------|
| `Validation/analytical.py` | Exact solution functions for all test cases |
| `Validation/convergence.py` | h-refinement and p-refinement study runner |
| `Validation/benchmarks.py` | De Vahl Davis and other benchmark runners |
| `Validation/views.py` | ConvergenceView, BenchmarkView, ValidationDashboard |
| `Validation/urls.py` | URL routes |
| `Validation/templates/convergence.html` | Convergence tables and plots |
| `Validation/templates/benchmark.html` | Benchmark comparison |
| `Validation/templates/dashboard.html` | Master validation dashboard |

---

## Definition of Done

- [x] h-refinement study for 1D and 2D shows rates ≈ 2.0 (L2) and 1.0 (H1) for P1 —
      `Validation.convergence.poisson_1d_convergence` /
      `poisson_2d_convergence`, measured rates 2.00 / 1.00 (1D) and 1.99 / 1.00 (2D)
      between the two finest meshes.
- [x] p-refinement study shows P1/P2/P3 converging at rates 2/3/4 (L2) —
      `Validation.convergence.p_refinement_study` on `n ∈ {4,8,16,32}`, measured
      L2 rates P1 ≈ 1.99, P2 ≈ 3.00, P3 ≈ 4.02.
- [x] Convection-diffusion validation confirms SUPG advantage at high Pe —
      *design refinement, see Completion Notes*: raw L2 error alone does not show
      a clean SUPG advantage (SUPG's added numerical diffusion can give it a
      slightly larger L2 error than an oscillatory Galerkin solution on the same
      mesh). The implemented check instead compares min/max nodal values: at
      `Pe=100`, `Pe_h=6.25 > 1`, Galerkin undershoots to `u_min ≈ -0.87` while
      SUPG stays bounded in `[0,1]` — the discrete-maximum-principle property
      that is SUPG's actual benefit.
- [x] De Vahl Davis benchmark page exists with Nusselt number comparison —
      `/Validation/benchmark/`, `Validation.benchmarks.de_vahl_davis_conduction`
      gives `Nu = 1.0000000000000080` for the Ra=0 limit (exact value 1.0).
- [x] Validation dashboard page links all studies — `/Validation/` (the
      pre-existing `'validation'` URL name, referenced by every page's navbar,
      now resolves here) shows 8 pass/fail checks covering every item above and
      links to `/Validation/convergence/` and `/Validation/benchmark/`.
- [x] All convergence plots use log-log scale with slope annotation —
      `_plot_series` in `Validation/convergence.py` renders `ax.loglog` with a
      dashed `O(h^p)` reference line and an "(observed ≈ X.XX)" rate annotation
      in the legend for every series.

## Completion Notes

### What was built

- **`Validation/` Django app** (no models, mirroring `TransientThermal/`):
  `analytical.py` (manufactured/exact solutions), `convergence.py` (h-refinement,
  p-refinement, and convection-diffusion convergence studies + log-log plotting),
  `benchmarks.py` (De Vahl Davis Ra=0 conduction limit + Nusselt number),
  `views.py` (`dashboard`, `convergence`, `benchmark`), `urls.py`, and templates
  `dashboard.html`, `convergence.html`, `benchmark.html`.
- **`_rate(prev_err, err, prev_h, h)`** implements the convergence-rate formula
  from `docs/science/governing_equations.md` §4
  (`rate = log(e_coarse/e_fine) / log(h_coarse/h_fine)`), used by every study.
- **h-refinement (4.2)**: `poisson_1d_convergence` and `poisson_2d_convergence`
  solve the manufactured-solution problems with P1 elements on
  `n ∈ {4,8,16,32,64}`, reporting L2/H1 errors and observed rates.
- **p-refinement (4.3)**: `p_refinement_study` calls `poisson_2d_convergence`
  for `orders=(1,2,3)`, confirming the `O(h^{p+1})` L2 rates.
- **Convection-diffusion convergence (4.4)**: `convection_diffusion_convergence`
  solves the 1D benchmark with both `'galerkin'` and `'supg'` for
  `Pe ∈ {0.1, 100}`, reporting L2 error and min/max nodal values.
- **De Vahl Davis (4.5)**: `de_vahl_davis_conduction(n)` solves the Ra=0
  pure-conduction limit (hot left wall, cold right wall, insulated top/bottom)
  and computes `Nu = -∫_{x=0} ∂u/∂x dy` via `assemble(grad(u_sol)[0]*ds(left))`.
  `REFERENCE_NUSSELT` tabulates the De Vahl Davis (1983) `Ra ∈ {0, 1e3, 1e4, 1e5,
  1e6}` reference values for context.
- **Dashboard (4.6)**: `/Validation/` re-runs every study at reduced resolution
  and renders 8 pass/fail checks against the expected rates/values; all 8 pass.
- **Architectural decision — repurposing the `'validation'` URL name**: the URL
  name `'validation'` was already referenced by every page's navbar
  (`{% url 'validation' %}`, 8+ templates) and previously pointed at a basic
  1D-only h-refinement check in `SteadyStateThermal`. That view, its URL, and
  its template were removed, and `Validation/urls.py` now defines
  `path('', dashboard, name='validation')` — so every existing navbar link now
  leads to the comprehensive new dashboard with **zero template edits**.
- **Documentation**: `docs/science/governing_equations.md` §4 gained an
  "Implementation Notes (`Validation/`)" subsection, and
  `docs/science/benchmark_references.md` gained a "Verified" note under each
  benchmark (sections 1–5) pointing to the corresponding `Validation`/
  `TransientThermal` function and its measured result.

### Test coverage

`tests/test_validation.py` — 38 new tests: analytical exact solutions, the
`_rate` helper, h-refinement and p-refinement structure/rates, convection-diffusion
convergence (low- and high-Pe behavior), plotting helpers, the De Vahl Davis
benchmark and Nusselt number, view-formatting helpers, and Django views
(dashboard, convergence, benchmark, and the repointed `'validation'` URL — including
a regression test that the old `/SteadyStateThermal/validation/` URL now 404s, and
that every existing page's navbar still resolves). Full suite: **117/117 passing**.

```text
117 passed in 17.56s
```

### Design decision — SUPG advantage measured via boundedness, not L2 error

Task 4.4 originally expected "SUPG should converge cleanly even at high Pe;
Galerkin will not on coarse meshes" to show up directly as a smaller L2 error
for SUPG. In practice, SUPG's artificial diffusion can make its L2 error
*comparable to or even larger than* an oscillatory Galerkin solution on the
same coarse mesh — L2 error alone doesn't capture the qualitative difference.
The dashboard check (and the convergence-page table) instead reports the
min/max nodal values: at `Pe=100` with `Pe_h > 1`, Galerkin produces large
unphysical undershoots (`u_min ≈ -0.87`, well outside the exact solution's
`[0,1]` range) while SUPG remains bounded — directly demonstrating the
discrete-maximum-principle property that motivates SUPG stabilization
(Goal 3).
