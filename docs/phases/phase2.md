# Phase 2 — Convection-Diffusion with SUPG Stabilization

**Duration:** 4 weeks  
**Status:** Not Started  
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

**This is the key figure** for your thesis chapter on stabilization.

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

- [ ] Convection-diffusion solver works for 1D and 2D meshes
- [ ] SUPG stabilization implemented and tested
- [ ] Peclet study page shows Galerkin oscillations vs SUPG correction side-by-side
- [ ] L2 error table comparing Galerkin vs SUPG vs exact solution
- [ ] Velocity field overlay shown on 2D solution plot
- [ ] Inflow/outflow BC handling works correctly
