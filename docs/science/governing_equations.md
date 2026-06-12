# Governing Equations Reference

This document is the mathematical reference for all PDE problems implemented or planned in FenicsWeb.

---

## 1. Steady-State Heat Conduction (Poisson Equation)

**Strong form:**
```
-∇·(k∇u) = f     in Ω
```

**Boundary conditions:**
- Dirichlet: `u = g`    on `Γ_D`
- Neumann:   `-k ∂u/∂n = q`    on `Γ_N`
- Robin:     `-k ∂u/∂n = h(u - u∞)`    on `Γ_R`

**Weak form:**
Find `u ∈ H¹(Ω)` such that for all `v ∈ H¹₀(Ω)`:
```
∫_Ω k∇u·∇v dx  =  ∫_Ω f·v dx  -  ∫_{Γ_N} q·v ds  -  ∫_{Γ_R} h(u-u∞)·v ds
```

In FEniCS:
```python
a = k * dot(grad(u), grad(v)) * dx + h * u * v * ds(robin_id)
L = f * v * dx - q * v * ds(neumann_id) + h * u_inf * v * ds(robin_id)
```

**Analytical 1D solution** (`-k u'' = f`, Dirichlet both ends):
```
u(x) = u_L + (u_R - u_L) * x/L  +  f/(2k) * x * (L - x)
```

---

## 2. Steady-State Convection-Diffusion

**Strong form:**
```
-k∇²u + b·∇u = f     in Ω
```

**Cell Peclet number:**
```
Pe_h = |b| * h / (2k)
```
- `Pe_h < 1`: diffusion-dominated — standard Galerkin works
- `Pe_h > 1`: convection-dominated — SUPG stabilization required

**Standard Galerkin weak form:**
```
∫_Ω k∇u·∇v dx  +  ∫_Ω (b·∇u) v dx  =  ∫_Ω f v dx
```

**SUPG stabilized weak form:**
Replace test function `v` with `v_SUPG = v + τ (b·∇v)`:
```
τ_K = h_K / (2|b|) * (coth(Pe_h) - 1/Pe_h)
```

**Analytical 1D solution** (`-u'' + Pe·u' = 0`, `u(0)=0`, `u(1)=1`):
```
u(x) = (exp(Pe·x) - 1) / (exp(Pe) - 1)
```

### Implementation Notes (`ConvectionDiffusion/solver.py`)

**Stabilization parameter.** The full SUPG formula `τ_K = h_K/(2|b|) · (coth(Pe_h) - 1/Pe_h)`
is singular as `Pe_h → 0` (i.e. `b → 0`, pure diffusion) because
`coth(Pe_h) - 1/Pe_h → ∞ - ∞`. The implementation uses the analytic small-`Pe_h`
limit `τ_K → h_K²/(12k)` (the first term of the Taylor expansion of
`coth(x) - 1/x = x/3 - x³/45 + …`) and switches between the two branches with
`ufl.conditional`:

```python
tau = conditional(
    gt(Pe, 1e-8),
    h / (2 * b_norm + DOLFIN_EPS) * (1 / tanh(Pe) - 1 / Pe),
    h * h / (12 * k),
)
```

Both terms of the bilinear form (diffusion and convection) are tested against
`v_SUPG = v + τ(b·∇v)`, matching the formulation in `phases/phase2.md`.

**Peclet number.** `compute_peclet_number(mesh, k, bx, by)` returns
`Pe_h = |b|·h_max/(2k)` using the largest cell diameter in the mesh, displayed
to the user on the solver page to indicate whether the problem is
diffusion-dominated (`Pe_h < 1`) or convection-dominated (`Pe_h > 1`).

**Inflow/outflow detection.** `mark_inflow_outflow(mesh, bx, by)` classifies
every exterior facet by the sign of `b·n` (velocity dotted with the outward
normal): `1` = inflow (`b·n < 0`, where a Dirichlet value is prescribed),
`2` = outflow (`b·n > 0`, natural/zero-flux), `3` = characteristic
(`b·n ≈ 0`, no-flux wall). This requires `mesh.init(dim-1, dim)` to populate
facet-normal connectivity before calling `Facet.normal()`.

**Peclet study benchmark.** `peclet_study_1d(pe_values, n)` solves the 1D
benchmark above for `Pe ∈ {1, 10, 50, 100}` with `k = 1/Pe`, `b = 1` on a
fixed 20-cell mesh, returning nodal values and L2 errors for both Galerkin
and SUPG against the exact solution. This is the figure used to verify SUPG
stabilization (`/ConvectionDiffusion/peclet/`).

---

## 3. Transient Heat Equation

**Strong form:**
```
ρ·cp · ∂u/∂t  =  ∇·(k∇u) + f     in Ω × [0,T]
```

With initial condition `u(x,0) = u₀(x)`.

**Backward Euler discretization** (1st-order, `θ=1`):
```
ρ·cp/dt * (u^n - u^(n-1)) = ∇·(k∇u^n) + f
```
Weak form: `a = ρcp/dt*(u,v) + k*(∇u,∇v)`,  `L = ρcp/dt*(u_n,v) + (f,v)`

**Crank-Nicolson discretization** (2nd-order, `θ=0.5`):
```
ρ·cp/dt * (u^n - u^(n-1)) = ∇·(k∇u_mid) + f,   u_mid = (u^n + u^(n-1))/2
```

**Temporal convergence rates:**
- Backward Euler: `O(dt¹)`
- Crank-Nicolson: `O(dt²)`

**Analytical solution** (homogeneous Dirichlet, no source, 1D):
```
u(x,t) = sin(πx) * exp(-π²·k·t / (ρ·cp))
```
(Use initial condition `u₀ = sin(πx)`)

### Implementation Notes (`TransientThermal/solver.py`)

**Time-stepping loop.** `solve_transient(...)` computes the number of steps
as `nsteps = max(1, round(T/dt))` and iterates `t = step*dt` for
`step = 1, …, nsteps`. An earlier `while t < T - 1e-10: t += dt` formulation
overshoots `T` whenever `T/dt` is not an integer (e.g. `dt=0.04, T=0.1` gives
`T/dt=2.5`, and the `while` loop ran 3 steps to `t=0.12`). The `round()`-based
step count guarantees `times[-1] == nsteps*dt` is the closest achievable
multiple of `dt` to `T`, with no overshoot.

**Weak forms.** Both schemes assemble `F` via UFL and split it with FEniCS's
`lhs`/`rhs` helpers:

```python
# Backward Euler
F = rho*cp/dt*(u - u_n)*v*dx + k*dot(grad(u), grad(v))*dx - Constant(f)*v*dx

# Crank-Nicolson
u_mid = 0.5*(u + u_n)
F = rho*cp/dt*(u - u_n)*v*dx + k*dot(grad(u_mid), grad(v))*dx - Constant(f)*v*dx

a, L = lhs(F), rhs(F)
```

**Initial conditions.** A fixed dictionary `INITIAL_CONDITIONS['1D'|'2D']`
maps the user-selectable keys `'sine'`, `'gaussian'`, `'constant'` to
hardcoded `Expression` strings. FEniCS `Expression` objects are JIT-compiled
as C++ at solve time, so the expression string is **never** built from raw
user input — only one of these three fixed strings is ever passed to
`Expression(...)`.

**Boundary conditions.** Homogeneous or non-homogeneous Dirichlet values are
applied on every domain edge (`left`/`right` in 1D; `+ bottom`/`top` in 2D)
using `near(x[i], bound)` matches against the mesh's own coordinate extents,
the same pattern used by `ConvectionDiffusion.solve_convection_diffusion`.

**Temporal convergence study.** `scheme_comparison_1d(dt_values, T, n=200)`
solves the `sin(πx)` benchmark (homogeneous Dirichlet, `f=0`) for each `dt`
with both schemes on a fixed `n=200` mesh (`h=0.005`) and compares the
solution at `t=T` against the exact decay solution via `errornorm(..., 'L2')`.
A mesh fine enough that the `O(h²)` spatial error (~2.5e-5) stays below the
temporal error is required to observe the theoretical rates; on a coarser
`n=100` mesh (`h=0.01`, spatial floor ~1e-4) the Crank-Nicolson rate degrades
from ≈2.0 toward ≈1.4 as `dt → 0`, because the spatial discretization error
dominates. With `n=200`, measured rates for `dt ∈ {0.02, 0.01, 0.005}` are
≈0.95/0.97 (Backward Euler, expected `O(dt¹)`) and ≈1.95/1.80
(Crank-Nicolson, expected `O(dt²)`). This study is reproduced at
`/TransientThermal/comparison/`.

**Visualization.** For 1D problems, `generate_transient_plotly(times,
snapshots, coords)` builds a single Plotly figure with one `go.Frame` per
time step plus a slider and play/pause buttons, so the full time history can
be scrubbed client-side without additional requests. `generate_transient_animation`
renders the same data as an animated GIF via `matplotlib.animation` (pillow
writer). For 2D problems, `generate_transient_snapshot_grid` plots the field
at evenly-spaced times as a row of `tripcolor` panels sharing one colorbar
(animation is skipped for 2D to avoid large in-page payloads).

**Resource guards.** `TransientThermal.forms.TransientForm` caps the mesh
resolution at `n ≤ 50` and validates `T/dt ≤ MAX_TIME_STEPS = 200` in
`clean()`, bounding solve time and the size of the generated Plotly/GIF
output for a web-exposed solver.

---

## 4. Finite Element Error Estimates

For a problem with solution `u ∈ H^{s+1}(Ω)` solved with P_p elements on mesh of size `h`:

| Norm | Rate |
|------|------|
| L2:   `||u - u_h||_0` | `O(h^{p+1})` |
| H1:   `||u - u_h||_1` | `O(h^p)` |
| Energy: `|u - u_h|_a` | `O(h^p)` |

**For P1 (linear):** L2 → O(h²), H1 → O(h¹)  
**For P2 (quadratic):** L2 → O(h³), H1 → O(h²)  
**For P3 (cubic):** L2 → O(h⁴), H1 → O(h³)

**Convergence rate computation:**
```
rate = log(e_coarse / e_fine) / log(h_coarse / h_fine)
```

### Implementation Notes (`Validation/`)

The `Validation` Django app is the verification framework required by Goal 6.
It is deliberately stateless (no models): every page recomputes its results
from scratch on each request, so every number on the site is reproducible by
re-running the corresponding function.

**`Validation/analytical.py`** provides the manufactured/exact solutions used
throughout: `exact_1d_poisson` (`u=x(1-x)`, Section 2 above),
`exact_2d_poisson` (`u=sin(πx)sin(πy)`, Section 3), and `source_2d_poisson`
(the matching right-hand side `f=2π²k·sin(πx)sin(πy)`).

**`Validation/convergence.py`** implements `_rate(prev_err, err, prev_h, h)`
exactly as the formula above, and three studies:

- `poisson_1d_convergence(mesh_sizes, k)` — h-refinement of `-k u''=2k` on
  `[0,1]` with P1 elements. Measured rates on `n ∈ {8,16,32,64}`: L2 → 2.00,
  H1 → 1.00 (matches theory to 4 decimal places).
- `poisson_2d_convergence(mesh_sizes, k, element_order)` — h-refinement of
  the 2D manufactured solution; `element_order` selects P1/P2/P3.
- `p_refinement_study(mesh_sizes, orders)` — calls `poisson_2d_convergence`
  for `orders=(1,2,3)`. Measured L2 rates on `n ∈ {4,8,16,32}`: P1 → 1.99,
  P2 → 3.00, P3 → 4.02.
- `convection_diffusion_convergence(pe_values, mesh_sizes)` — solves the 1D
  convection-diffusion benchmark (Section 4, ConvectionDiffusion solver)
  with both `'galerkin'` and `'supg'` and reports L2 error **and** the
  min/max nodal value. The min/max pair is the more informative diagnostic
  at high Pe: at `Pe=100`, `Pe_h=6.25` (`n=8`), standard Galerkin undershoots
  to `u_min ≈ -0.87` (a 87% spurious oscillation below the physical lower
  bound `u=0`), while SUPG remains within `[0,1]` at every mesh size tested.
  Raw L2 error is *not* a reliable indicator of SUPG's benefit here — SUPG's
  added numerical diffusion can give it a marginally larger L2 error than an
  oscillatory Galerkin solution on the same mesh, even though the Galerkin
  solution is qualitatively wrong (violates the discrete maximum principle).

All three studies feed `generate_h_refinement_plot`,
`generate_p_refinement_plot`, and `generate_convection_diffusion_plot`,
which render log-log error-vs-`h` plots with a dashed `O(h^p)` reference
line and the *observed* rate (computed via `_rate` between the coarsest and
finest mesh) annotated directly in the legend.

**`Validation/benchmarks.py`** implements the De Vahl Davis cavity
(Section 1) in its Ra=0 (pure-conduction) limit:
`de_vahl_davis_conduction(n)` solves `-∇²u=0` on `[0,1]²` with `u=1` on the
left wall, `u=0` on the right wall, and natural (zero-flux) conditions on
top/bottom — i.e. those boundaries are simply omitted from the weak form.
The average Nusselt number is computed as
`Nu = -∫_{x=0} ∂u/∂x dy = -assemble(grad(u_sol)[0]*ds(left))`, using a
`MeshFunction` to mark the left-wall facets. The measured value is
`Nu = 1.0000000000000080` for `n=32` — i.e. exact to 14 significant digits,
since the discrete solution of `-∇²u=0` with these boundary conditions is
the linear function `u=1-x`, which P1 elements represent exactly.

**`Validation/views.py`** provides three pages: `/Validation/convergence/`
(all four studies above, full mesh-size ranges), `/Validation/benchmark/`
(De Vahl Davis), and `/Validation/` — the validation dashboard, which
re-runs each study at reduced resolution (for response time) and renders a
table of pass/fail checks against the expected rates/values documented here
and in `docs/science/benchmark_references.md`.
