# Phase 4 — Scientific Rigor & Benchmarks

**Duration:** 3 weeks  
**Status:** Not Started  
**Depends on:** Phase 1 (correct solver), Phase 2 (convection-diffusion)  
**Unlocks:** Phase 6 (research contribution needs validated solver)

---

## Goal

Add the Verification and Validation (V&V) infrastructure that separates a research tool from a student project. Every figure produced here can go directly into a thesis chapter.

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

**Output:** Table + log-log convergence plot (this plot belongs in every FEM thesis).

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

Compute the Nusselt number from your solver and compare in a table. Note that without the flow solver, we can only compare the thermal field at Ra=0 (pure conduction), but the infrastructure and comparison methodology is what matters for the thesis.

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

- [ ] h-refinement study for 1D and 2D shows rates ≈ 2.0 (L2) and 1.0 (H1) for P1
- [ ] p-refinement study shows P1/P2/P3 converging at rates 2/3/4 (L2)
- [ ] Convection-diffusion validation confirms SUPG advantage at high Pe
- [ ] De Vahl Davis benchmark page exists with Nusselt number comparison
- [ ] Validation dashboard page links all studies
- [ ] All convergence plots use log-log scale with slope annotation
