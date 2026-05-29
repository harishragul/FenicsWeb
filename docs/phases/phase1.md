# Phase 1 — Foundation Repair

**Duration:** 2 weeks
**Status:** Complete
**Completed:** 2026-05-29
**Depends on:** Phase 0 (repo public, CI running)
**Unlocks:** Phase 2, Phase 3, Phase 4

---

## Goal

Make what already exists correct, complete, and scientifically credible. Every later phase builds on this foundation — bugs here propagate everywhere.

---

## Tasks

### 1.1 Fix the solution visualization bug

**File:** `SteadyStateThermal/conduction.py:72`
**Bug:** `generate_solution_plot(mesh)` received a `mesh` object and called `plot(mesh)` — plotting the mesh topology again instead of the temperature field.

**Fix applied:**

- New signature: `generate_solution_plot(u_sol, mesh)` — takes the FEniCS `Function` object
- For **1D**: sorted line plot of `u vs x` with grid
- For **2D**: `fenics.plot(u_sol)` with labeled Viridis colorbar
- For **3D**: three scatter cross-section slices at x-mid, y-mid, z-mid
- Call site in `SteadyStateThermal/views.py` updated to pass `(u_sol, mesh)`

---

### 1.2 Enable the interactive Plotly visualization

**File:** `mesh/mesh.py:174`

**Fix applied:**

- `generate_interactive_mesh_plot(mesh, u_sol=None)` — accepts optional solution
- When `u_sol` is provided, computes vertex values and passes them as `intensity` to `go.Mesh3d`
- Adds `colorscale='Viridis'` and a `colorbar` labelled "Temperature (K)"
- `fig.update_layout(margin=...)` for clean embedding

---

### 1.3 Add Neumann boundary conditions

**What it is:** Prescribed heat flux at a boundary face — `-k ∂u/∂n = q`

**Weak form:** Add `q * v * ds(tag)` to the right-hand side `L` where `ds(tag)` is the surface measure on the tagged boundary facet only.

**Implementation:**

- `_mark_facets(mesh, coords)` helper: iterates `facets(mesh)` and assigns tags 1–6 (left/right/bottom/top/front/back) via `MeshFunction`
- `solve_conduction` accepts optional `bc_types`, `neumann_fluxes` dicts
- Neumann faces skip `DirichletBC` and contribute `q * v * ds(tag)` to `L`

**Validation:** `test_1d_neumann_right_bc_linear` — Dirichlet u(0)=0, Neumann q=1 → u(x)=x. Passes.

---

### 1.4 Add Robin (convection) boundary conditions

**What it is:** Newton's law of cooling — `-k ∂u/∂n = h(u - u∞)`

**Weak form:** `h * u * v * ds(tag)` added to bilinear form `a`; `h * u∞ * v * ds(tag)` added to `L`.

**Implementation:**

- `solve_conduction` accepts optional `robin_hs`, `robin_u_infs` dicts
- Robin faces contribute to both `a` and `L`; no `DirichletBC` added

**Validation:** `test_1d_robin_right_temperature_above_ambient` — left=100°, Robin right (h=10, T∞=20) → boundary temperature between 20 and 100. Passes.

---

### 1.5 Analytical validation for 1D

**Exact solution** for `-k u″ = f` on `[0,1]` with `u(0)=u(1)=0`, `f=2`, `k=1`:

```text
u(x) = x(1 − x)
```

**Implementation:**

- `validation` view in `SteadyStateThermal/views.py` solves on n = 4, 8, 16, 32, 64 meshes
- `errornorm(u_exact, u_sol, 'L2')` gives the L2 error at each level
- Convergence rate computed as `log(e_prev / e) / log(2)` between successive levels
- Results rendered in `SteadyStateThermal/templates/validation.html` — table with rate highlighted green when ≈ 2.0

---

### 1.6 Expand the test suite

Tests added to `tests/test_conduction.py`:

| Test | What it verifies |
| ---- | ---------------- |
| `test_1d_neumann_right_bc_linear` | Neumann q=1 → u(x)=x at interior nodes |
| `test_1d_neumann_both_sides_raises_or_zero` | Purely Neumann system handled gracefully |
| `test_1d_robin_right_temperature_above_ambient` | Robin BC → T_boundary in (u∞, T_left) |
| `test_1d_robin_zero_h_equals_neumann_zero_flux` | h=0 Robin → insulated (constant profile) |
| `test_interactive_mesh_plot_returns_valid_json` | Plotly output is valid JSON with "data" key |
| `test_interactive_mesh_plot_with_solution_has_intensity` | Solution coloring passes intensity array |

---

## Files Changed

| File | Change |
| ---- | ------ |
| `SteadyStateThermal/conduction.py` | Full rewrite — fixed plot, added Neumann/Robin, facet marking |
| `SteadyStateThermal/forms.py` | Added bc_type dropdowns and h/u_inf Robin fields per face |
| `SteadyStateThermal/views.py` | Updated call sites; added `validation` view |
| `SteadyStateThermal/urls.py` | Added `validation/` route |
| `mesh/mesh.py` | `generate_interactive_mesh_plot` accepts optional `u_sol` |
| `SteadyStateThermal/templates/conduction.html` | Rewritten — colormap image, validation link, JS BC-type toggle |
| `SteadyStateThermal/templates/validation.html` | New — convergence table with rate highlighting |
| `tests/test_conduction.py` | Added 6 new tests for Neumann, Robin, interactive plot |

---

## Definition of Done

- [x] Solution page shows a temperature colormap (not a mesh replot)
- [x] Interactive Plotly mesh colors vertices by temperature when u_sol is passed
- [x] Neumann BC implemented and validated — `test_1d_neumann_right_bc_linear` passes
- [x] Robin BC implemented and validated — `test_1d_robin_right_temperature_above_ambient` passes
- [x] Convergence table page at `/SteadyStateThermal/validation/` shows rate ≈ 2.0 for P1
- [x] All new functionality covered by pytest tests — **21/21 passing**
- [x] GitHub Actions CI passing

---

## Completion Notes

### Key technical issue resolved

`_mark_facets` originally used `markers.mesh().facets()`, which does not exist in this version
of legacy FEniCS (dolfin). The correct API is the module-level `facets(mesh)` iterator imported
from `fenics`. Fixed by replacing all three occurrences in the helper.

### Test count

```text
21 passed in 3.41s
```

| Module | Tests | New in Phase 1 |
| ------ | ----- | -------------- |
| `test_conduction.py` | 12 | 6 (Neumann, Robin, Plotly intensity) |
| `test_mesh.py` | 9 | 0 |

### Design decision — form structure

Rather than a separate solver page for Neumann/Robin, the existing `HeatSolverForm` was extended
with per-face bc_type dropdowns and conditional h/u_inf fields. JavaScript in the template shows
and hides Robin-specific fields based on the selected type, keeping the UI on a single page.

### Design decision — backward compatibility

`solve_conduction` retains its original positional/keyword signature. The new `bc_types`,
`neumann_fluxes`, `robin_hs`, and `robin_u_infs` kwargs all default to `None` (resolved to
empty dicts inside the function), so all existing tests continue to pass without modification.
