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

This is a clean temporal convergence result for the thesis.

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

- [ ] Backward Euler time-stepping works for 1D and 2D
- [ ] Crank-Nicolson time-stepping works for 1D and 2D
- [ ] Snapshots saved to database
- [ ] Browser time scrubber re-renders plot via AJAX without page reload
- [ ] Animation GIF/MP4 export works
- [ ] Temporal convergence table shows BE ≈ 1st-order, CN ≈ 2nd-order
- [ ] Scheme comparison page built with exact-solution overlay
