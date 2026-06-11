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
