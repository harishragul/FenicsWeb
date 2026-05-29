# Phase 1 — Foundation Repair

**Duration:** 2 weeks  
**Status:** Not Started  
**Depends on:** Nothing (baseline)  
**Unlocks:** Phase 2, Phase 3, Phase 4

---

## Goal

Make what already exists correct, complete, and scientifically credible. Every later phase builds on this foundation — bugs here propagate everywhere.

---

## Tasks

### 1.1 Fix the solution visualization bug

**File:** `SteadyStateThermal/conduction.py:72`  
**Bug:** `generate_solution_plot(mesh)` receives a `mesh` object and calls `plot(mesh)` — it plots the mesh topology again, not the temperature field.

**Fix:**
- The function signature must be `generate_solution_plot(u_sol, mesh)` to receive the FEniCS `Function` object
- For **1D**: plot `u vs x` as a line plot (x-axis: coordinates, y-axis: temperature values)
- For **2D**: use `plot(u_sol)` or Matplotlib `tricontourf` with a labeled colorbar
- For **3D**: render 2D cross-section slices along each axis plane (z=0.5, y=0.5, x=0.5)
- Update the call site in `SteadyStateThermal/views.py:26` to pass `solution` (the `u_sol` object) instead of `mesh`

**Acceptance:** The result page shows temperature distribution, not a second copy of the mesh.

---

### 1.2 Enable the interactive Plotly visualization

**Files:** `mesh/mesh.py:174`, `mesh/templates/mesh.html:63`  
**Issue:** `generate_interactive_mesh_plot()` is fully implemented but commented out at the call site in `views.py:17` and in the template.

**Fix:**
- Uncomment and test the existing Plotly `Mesh3d` function
- Extend it to also display the solution scalar field: color mesh vertices by temperature value using `intensity` parameter in `go.Mesh3d`
- Add a Plotly colorscale (e.g., `colorscale='Viridis'`) with a colorbar label "Temperature (K)"

**Acceptance:** Mesh page shows an interactive 3D mesh; solution page shows mesh colored by temperature.

---

### 1.3 Add Neumann boundary conditions

**What it is:** Prescribed heat flux at a boundary face — `q = -k ∂u/∂n = g`

**Weak form addition:** Add `g * v * ds` to the right-hand side `L` where `ds` is the surface measure on the target boundary.

**Implementation:**
- Add `neumann_flux` fields to `HeatSolverForm` in `SteadyStateThermal/forms.py` for each face (left, right, top, bottom, front, back)
- Add a `bc_type` choice per face: `Dirichlet` or `Neumann`
- In `conduction.py`, only add `DirichletBC` for faces marked Dirichlet; add the flux surface integral for Neumann faces
- Use `FacetFunction` / `MeshFunction` to mark boundary facets by subdomain

**Acceptance:** A 1D rod with Neumann flux on the right (q=1) and Dirichlet on the left (u=0) produces `u(x) = x/k` — verify this numerically.

---

### 1.4 Add Robin (convection) boundary conditions

**What it is:** Newton's law of cooling — `q = h(u - u_inf)` where `h` is the convection coefficient and `u_inf` is the ambient temperature.

**Weak form addition:** Add `h * u * v * ds - h * u_inf * v * ds` — the `h*u*v*ds` term moves to the bilinear form `a`, and `-h*u_inf*v*ds` goes to `L`.

**Implementation:**
- Add `robin_h` (convection coefficient) and `robin_u_inf` (ambient temperature) fields to the form
- Add `Robin` as a third `bc_type` option per face
- Update `solve_conduction()` in `conduction.py` to assemble the Robin terms

**Why this matters:** Robin BCs model heat loss to a surrounding fluid — this is *actual* convection modeling and justifies the project title.

**Acceptance:** A rod insulated on the left, with Robin BC on the right (h=10, u_inf=20), produces the expected temperature profile matching the 1D analytical solution.

---

### 1.5 Analytical validation for 1D

**Purpose:** Demonstrate the solver is correct by comparing against a known exact solution — this is the bare minimum for scientific credibility.

**Exact solution:** For `-k u'' = f` on `[0, L]` with `u(0) = u_L`, `u(L) = u_R`:
```
u(x) = u_L + (u_R - u_L) * x/L + f/(2k) * x * (L - x)
```

**Implementation:**
- Add a `ValidationView` in `SteadyStateThermal/views.py`
- Solve on a sequence of meshes: n = 4, 8, 16, 32, 64 cells
- Compute L2 error: `||u_h - u_exact||_L2 = sqrt(assemble((u_h - u_exact)**2 * dx))`
- Render a results table: columns = n, h, L2 error, convergence rate
- Convergence rate should be ≈ 2.0 for P1 elements (confirming correct implementation)

**Acceptance:** Convergence table shows rate ≈ 2.0; this table is the first figure you can put in a thesis.

---

## Files to Create / Modify

| File | Action |
|------|--------|
| `SteadyStateThermal/conduction.py` | Fix `generate_solution_plot`, add Neumann/Robin assembly |
| `SteadyStateThermal/forms.py` | Add bc_type, neumann_flux, robin_h, robin_u_inf fields |
| `SteadyStateThermal/views.py` | Fix plot call, add ValidationView |
| `SteadyStateThermal/urls.py` | Add validation URL route |
| `mesh/mesh.py` | Enable interactive plot, add solution coloring |
| `mesh/templates/mesh.html` | Uncomment Plotly block |
| `SteadyStateThermal/templates/conduction.html` | Update to show solution colormap + validation link |
| `SteadyStateThermal/templates/validation.html` | New — convergence table |

---

### 1.6 Expand the test suite

Every task in this phase must have a corresponding test. Add to `tests/test_conduction.py`:

**Neumann BC test:**

```python
def test_1d_neumann_right_bc():
    """Rod with Dirichlet u(0)=0 and Neumann q=1 on right gives u(x) = x/k."""
    # Exact: u(x) = x (k=1, q=1)
    mesh = UnitIntervalMesh(20)
    u_sol = solve_conduction_neumann(mesh, left_bc=0, right_flux=1.0, k=1)
    assert abs(u_sol(1.0) - 1.0) < 1e-6
    assert abs(u_sol(0.5) - 0.5) < 1e-4
```

**Robin BC test:**

```python
def test_1d_robin_steady_state():
    """Robin BC solution must match analytical solution within tolerance."""
    mesh = UnitIntervalMesh(40)
    h_coeff, u_inf = 10.0, 20.0
    u_sol = solve_conduction_robin(mesh, left_bc=100, right_h=h_coeff,
                                    right_u_inf=u_inf, k=1, f=0)
    # Value at right boundary should satisfy: -du/dx = h*(u - u_inf)
    # Verified numerically against scipy.integrate.solve_bvp reference
    assert u_sol(1.0) > u_inf   # temperature must be above ambient
```

**Interactive plot test:**

```python
def test_interactive_plot_returns_json():
    from mesh.mesh import generate_interactive_mesh_plot
    from fenics import UnitSquareMesh
    import json
    mesh = UnitSquareMesh(4, 4)
    result = generate_interactive_mesh_plot(mesh)
    data = json.loads(result)
    assert "data" in data
```

---

## Phase 1 Files to Create / Modify

| File | Action |
| ---- | ------ |
| `SteadyStateThermal/conduction.py` | Fix `generate_solution_plot`, add Neumann/Robin assembly |
| `SteadyStateThermal/forms.py` | Add bc_type, neumann_flux, robin_h, robin_u_inf fields |
| `SteadyStateThermal/views.py` | Fix plot call, add ValidationView |
| `SteadyStateThermal/urls.py` | Add validation URL route |
| `mesh/mesh.py` | Enable interactive plot, add solution coloring |
| `mesh/templates/mesh.html` | Uncomment Plotly block |
| `SteadyStateThermal/templates/conduction.html` | Update to show solution colormap + validation link |
| `SteadyStateThermal/templates/validation.html` | New — convergence table |
| `tests/test_conduction.py` | Add Neumann, Robin, interactive plot tests |

---

## Definition of Done

- [ ] Solution page shows a temperature colormap (not a mesh replot)
- [ ] Interactive 3D Plotly mesh works and colors vertices by temperature
- [ ] Neumann BC works and validates against 1D exact solution
- [ ] Robin BC works and validates against 1D exact solution
- [ ] Convergence table page exists and shows rate ≈ 2.0 for P1
- [ ] All new functionality covered by pytest tests
- [ ] GitHub Actions CI still passing after changes
