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

## 4. Inverse Heat Conduction

**Forward problem:** Given `k(x)`, find `u`:
```
-∇·(k∇u) = f    in Ω
```

**Inverse problem:** Given noisy measurements `{u(xᵢ)}_obs`, find `k(x)` minimizing:
```
J(k) = Σᵢ (u(xᵢ) - u_obs_i)²  +  α ∫_Ω |∇k|² dx
```

The second term is Tikhonov regularization (H¹ semi-norm of k).

**Adjoint equation** (for gradient computation):
```
-∇·(k∇p) = -2 Σᵢ (u(xᵢ) - u_obs_i) δ(x - xᵢ)
```
where `p` is the adjoint variable (Lagrange multiplier).

**Gradient of J with respect to k:**
```
dJ/dk = ∇u · ∇p  +  2α (-∇²k)
```

This is computed automatically by `dolfin-adjoint`.

---

## 5. Finite Element Error Estimates

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
