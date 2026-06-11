# Phase 3 — Transient Heat Equation

**Duration:** 4 weeks  
**Status:** Not Started  
**Depends on:** Phase 1 (correct BCs and plots)  
**Unlocks:** Phase 6C (UQ on transient problems)

---

## Goal

Add time-dependency to the heat equation. This is essential for modeling real thermal processes (heating up, cooling down, thermal shock). The scientific contribution is a comparison of two time-stepping schemes (Backward Euler vs Crank-Nicolson) and their accuracy/stability tradeoffs.

---

## Governing Equation

```
ρ·cp · ∂u/∂t  =  ∇·(k∇u) + f     in Ω × [0, T]
```

With:
- Initial condition: `u(x, 0) = u_0(x)`
- Boundary conditions: same as steady-state (Dirichlet, Neumann, Robin)

Parameters:
- `ρ` = density (kg/m³)
- `cp` = specific heat capacity (J/kg·K)
- `k` = thermal conductivity (W/m·K)
- `T` = end time (s)
- `dt` = time step (s)

---

## Time Discretization Schemes

### Scheme 1: Backward Euler (1st-order)

```
ρ·cp · (u^n - u^(n-1)) / dt  =  ∇·(k∇u^n) + f
```

Weak form (fully implicit — use `u` = current unknown, `u_n` = previous):
```python
F = rho*cp/dt * (u - u_n) * v * dx + k * dot(grad(u), grad(v)) * dx - f * v * dx
```

**Stability:** Unconditionally stable (any `dt`), but only 1st-order accurate in time.

### Scheme 2: Crank-Nicolson (2nd-order) — Primary Scheme

```
ρ·cp · (u^n - u^(n-1)) / dt  =  ∇·(k·∇((u^n + u^(n-1))/2)) + f
```

Weak form:
```python
u_mid = 0.5 * (u + u_n)
F = rho*cp/dt * (u - u_n) * v * dx + k * dot(grad(u_mid), grad(v)) * dx - f * v * dx
```

**Stability:** Unconditionally stable, 2nd-order accurate in time — preferred scheme.

---

## Tasks

### 3.1 Create the `TransientThermal` Django app

```
python manage.py startapp TransientThermal
```

```
TransientThermal/
├── __init__.py
├── apps.py
├── forms.py          # TransientForm: rho, cp, k, f, T, dt, scheme choice, u0
├── solver.py         # BackwardEuler, CrankNicolson solvers
├── models.py         # TimeSeriesResult: job_id, timestep, solution_data (JSON)
├── views.py          # TransientView, AnimationView, SchemeComparisonView
├── urls.py
└── templates/
    ├── transient.html
    ├── animation.html
    └── scheme_comparison.html
```

---

### 3.2 Implement the time-stepping loop

**File:** `TransientThermal/solver.py`

```python
def solve_transient(mesh, rho, cp, k, f, T, dt, scheme, u0_expr, bcs):
    V = FunctionSpace(mesh, "P", 1)
    u = TrialFunction(V)
    v = TestFunction(V)
    u_n = interpolate(Expression(u0_expr, degree=1), V)

    snapshots = [(0.0, u_n.compute_vertex_values(mesh).tolist())]
    t = dt

    while t <= T + 1e-10:
        if scheme == "backward_euler":
            F = rho*cp/dt*(u - u_n)*v*dx + k*dot(grad(u), grad(v))*dx - Constant(f)*v*dx
        else:  # crank_nicolson
            u_mid = 0.5*(u + u_n)
            F = rho*cp/dt*(u - u_n)*v*dx + k*dot(grad(u_mid), grad(v))*dx - Constant(f)*v*dx

        a, L = lhs(F), rhs(F)
        u_sol = Function(V)
        solve(a == L, u_sol, bcs)

        snapshots.append((round(t, 10), u_sol.compute_vertex_values(mesh).tolist()))
        u_n.assign(u_sol)
        t += dt

    return snapshots
```

---

### 3.3 Store time series in the database

**File:** `TransientThermal/models.py`

```python
class TransientJob(models.Model):
    name = models.CharField(max_length=200)
    created_at = models.DateTimeField(auto_now_add=True)
    scheme = models.CharField(max_length=50)
    dt = models.FloatField()
    T = models.FloatField()
    mesh_str = models.TextField()

class TransientSnapshot(models.Model):
    job = models.ForeignKey(TransientJob, on_delete=models.CASCADE)
    time = models.FloatField()
    solution_json = models.TextField()   # JSON-encoded vertex values
```

---

### 3.4 Browser-based time scrubber

**Template:** `TransientThermal/templates/transient.html`

Use a JavaScript range slider (`<input type="range">`) + AJAX to fetch the snapshot at a given timestep and re-render the plot without a page reload:

```javascript
slider.addEventListener('input', function() {
    const t = parseFloat(this.value);
    fetch(`/transient/snapshot/?job_id=${jobId}&t=${t}`)
        .then(r => r.json())
        .then(data => updatePlot(data.solution, data.coordinates));
});
```

Add a `SnapshotAPIView` in `views.py` that returns the snapshot JSON for a given `(job_id, t)`.

---

### 3.5 Animation export

Generate an animated GIF or MP4 of the temperature field over time:

```python
import matplotlib.animation as animation

def generate_animation(snapshots, coordinates):
    fig, ax = plt.subplots()
    ...
    def update(frame):
        ax.clear()
        # tricontourf plot for frame
    anim = animation.FuncAnimation(fig, update, frames=len(snapshots))
    anim.save('static/transient_animation.gif', writer='pillow', fps=10)
```

Add a download link on the results page.

---

### 3.6 Scheme comparison page

**Purpose:** Side-by-side accuracy comparison of Backward Euler vs Crank-Nicolson.

**Test problem:** 1D rod, initial condition `u(x,0) = sin(πx)`, no source, Dirichlet `u=0` at both ends.

**Exact solution:** `u(x,t) = sin(πx) · exp(-π²·k·t / (ρ·cp))`

**Show:**
- Temperature at `x = 0.5` vs time for both schemes and the exact solution
- L2 error at final time `T` for both schemes
- Convergence in `dt`: halve the time step, error should halve (BE) or quarter (CN)

This is a clean temporal convergence result demonstrating the expected first- and second-order accuracy of Backward Euler and Crank-Nicolson (Goal 4).

---

### 3.7 Analytical validation for transient

Using the same test problem as 3.6, build an automated validation table:

| Scheme | dt | L2 error at T | Rate |
|--------|----|---------------|------|
| Backward Euler | 0.1 | ... | — |
| Backward Euler | 0.05 | ... | ≈ 1.0 |
| Crank-Nicolson | 0.1 | ... | — |
| Crank-Nicolson | 0.05 | ... | ≈ 2.0 |

Expected rates: 1.0 for BE, 2.0 for CN — this verifies both implementations.

---

## Files to Create

| File | Purpose |
|------|---------|
| `TransientThermal/__init__.py` | App init |
| `TransientThermal/apps.py` | App config |
| `TransientThermal/forms.py` | Transient problem parameters form |
| `TransientThermal/solver.py` | Backward Euler + Crank-Nicolson solvers |
| `TransientThermal/models.py` | TransientJob, TransientSnapshot models |
| `TransientThermal/views.py` | TransientView, AnimationView, SnapshotAPIView, SchemeComparisonView |
| `TransientThermal/urls.py` | URL routes |
| `TransientThermal/templates/transient.html` | Solve page with time slider |
| `TransientThermal/templates/animation.html` | Animation download page |
| `TransientThermal/templates/scheme_comparison.html` | BE vs CN comparison |

**Files to modify:**
| File | Change |
|------|--------|
| `FenicsWeb/urls.py` | Include TransientThermal URLs |
| `FenicsWeb/settings.py` | Add TransientThermal to INSTALLED_APPS |

---

## Definition of Done

- [x] Backward Euler time-stepping works for 1D and 2D
- [x] Crank-Nicolson time-stepping works for 1D and 2D
- [x] Snapshots saved to database — *design change, see Completion Notes*: snapshots are
      embedded directly in the response (Plotly frames for 1D, a static snapshot grid for
      2D) instead of a `TransientJob`/`TransientSnapshot` DB schema, per CLAUDE.md's
      "no unnecessary persistence/backend complexity" guidance.
- [x] Browser time scrubber re-renders plot via AJAX without page reload — *design
      change*: 1D uses a client-side Plotly slider with all frames embedded in the page
      (no AJAX/DB round trip needed); 2D uses a static multi-panel snapshot grid instead
      of a scrubber (animated 2D payloads were judged too large for a web response).
- [x] Animation GIF/MP4 export works — implemented for 1D (`generate_transient_animation`,
      pillow GIF writer). Skipped for 2D, consistent with the snapshot-grid decision above.
- [x] Temporal convergence table shows BE ≈ 1st-order, CN ≈ 2nd-order
- [x] Scheme comparison page built with exact-solution overlay

## Completion Notes

### What was built

- **`TransientThermal/` Django app** mirroring `ConvectionDiffusion/`: `forms.py`,
  `solver.py`, `views.py`, `urls.py`, `templates/transient.html`,
  `templates/scheme_comparison.html`. No `models.py` content (see DoD notes above).
- **`solve_transient(dimension, mesh, rho, cp, k, f, T, dt, scheme, initial_condition,
  left_bc, right_bc, bottom_bc, top_bc)`** — single solver function for 1D
  (`UnitIntervalMesh`) and 2D (`UnitSquareMesh`) problems, returning
  `(times, snapshots, coords, u_final)`. Both Backward Euler and Crank-Nicolson weak forms
  are assembled with `lhs`/`rhs` exactly as drafted in §3.2, with Dirichlet BCs applied on
  every domain edge using the mesh's own coordinate bounds (same pattern as
  `ConvectionDiffusion.solve_convection_diffusion`).
- **No-overshoot time-stepping**: `nsteps = max(1, round(T/dt))`, iterating
  `t = step*dt` for `step = 1..nsteps`. The `while t <= T + 1e-10` loop drafted in §3.2
  overshoots `T` whenever `T/dt` is non-integer (e.g. `dt=0.04, T=0.1` → 3 steps to
  `t=0.12` instead of 2 steps to `t=0.08`), which corrupted the convergence comparison
  during development — fixed before any UI/test code was written.
- **Initial conditions**: a fixed `INITIAL_CONDITIONS['1D'|'2D']` dict maps the
  user-selectable keys `'sine'`, `'gaussian'`, `'constant'` to hardcoded `Expression`
  strings — never built from raw user input (FEniCS `Expression` is JIT-compiled C++).
- **`exact_transient_1d(x, t, k, rho, cp)`** — analytical decay solution
  `sin(πx)·exp(-π²kt/(ρcp))` for the 1D rod benchmark (§3.6/3.7).
- **`scheme_comparison_1d(dt_values, T, k, rho, cp, n)`** — the temporal-convergence
  benchmark. On a fixed `n=200` mesh (`h=0.005`, chosen so the `O(h²)` spatial error stays
  below the temporal error), measured rates for `dt ∈ {0.02, 0.01, 0.005}` are
  ≈0.95/0.97 for Backward Euler (expected `O(dt¹)`) and ≈1.95/1.80 for Crank-Nicolson
  (expected `O(dt²)`).
- **Visualization**: `generate_transient_plotly` (1D Plotly figure with one `go.Frame` per
  time step + slider + play/pause buttons), `generate_transient_animation` (1D animated
  GIF via `matplotlib.animation` + pillow), `generate_transient_snapshot_grid` (2D
  `tripcolor` panel grid at evenly-spaced times), `generate_scheme_comparison_plot`
  (two-panel BE/CN/exact comparison + log-log convergence plot with `O(dt)`/`O(dt²)`
  reference slopes).
- **Resource guards**: `TransientForm` caps mesh resolution at `n ≤ 50` and validates
  `T/dt ≤ MAX_TIME_STEPS = 200` in `clean()`.
- **UI integration**: navbar links ("Transient", "Scheme Comparison") added across all
  six existing pages; new `/TransientThermal/` solve page and
  `/TransientThermal/comparison/` benchmark page.
- **Documentation**: `docs/science/governing_equations.md` §3 gained an "Implementation
  Notes" subsection covering the time-stepping fix, weak forms, IC safety, and the
  convergence-study mesh-resolution requirement.

### Test coverage

`tests/test_transient.py` — 31 new tests (solver structure, no-overshoot time-stepping,
steady-state preservation for both schemes in 1D/2D, initial conditions, BE/CN accuracy
vs. the exact decay solution, temporal convergence rates, plotting/export helpers, form
validation, Django views). Full suite: **79/79 passing**.

```text
79 passed in 4.09s
```

### Design decision — DB-free snapshots and 2D visualization

Per CLAUDE.md, the project must remain a compact research platform rather than acquire
unnecessary backend/database complexity. Tasks 3.3/3.4 (`TransientJob`/`TransientSnapshot`
models + AJAX `SnapshotAPIView`) were replaced with embedding all time-step data directly
in the page response: 1D snapshots become Plotly `frames`/`sliders` (client-side scrubbing,
no extra requests), and 2D snapshots become a static `tripcolor` panel grid (animated 2D
payloads were judged too large to embed). This preserves every numerical-method and
verification requirement of Goal 4 while avoiding new models, migrations, or API
endpoints.
