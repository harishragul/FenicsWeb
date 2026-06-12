# Benchmark Problems and Reference Values

---

## 1. De Vahl Davis Natural Convection Cavity (1983)

**Reference:** De Vahl Davis, G. (1983). Natural convection of air in a square cavity: a benchmark numerical solution. *International Journal for Numerical Methods in Fluids*, 3(3), 249–264.

**Problem:** Square cavity `[0,1]²`, differentially heated vertical walls.
- Left wall: `T = 1` (hot)
- Right wall: `T = 0` (cold)
- Top and bottom: adiabatic (`∂T/∂n = 0`)

**Reference Nusselt numbers** (average over hot wall):

| Ra | Nu (De Vahl Davis) | Nu (Markatos & Pericleous 1984) |
|----|--------------------|----------------------------------|
| 10³ | 1.118 | 1.108 |
| 10⁴ | 2.243 | 2.201 |
| 10⁵ | 4.519 | 4.430 |
| 10⁶ | 8.800 | 8.754 |

**Note for this project:** FenicsWeb solves only the thermal diffusion problem (Ra=0 limit — no buoyancy-driven flow). The thermal field comparison at Ra=0 (pure conduction) gives `Nu = 1.0`. The full Boussinesq flow would require coupling the Navier-Stokes equations — a future extension.

**Verified:** `Validation.benchmarks.de_vahl_davis_conduction(n)` (`/Validation/benchmark/`) solves this Ra=0 limit on a 32×32 mesh and computes `Nu = -∫_{x=0} ∂u/∂x dy` via `assemble(grad(u_sol)[0]*ds(left))`. Measured `Nu = 1.0000000000000080`, exact to 14 significant digits.

---

## 2. 1D Steady Conduction — Analytical Benchmark

**Problem:** `-k u'' = f` on `[0, L]`, `u(0) = u_L`, `u(L) = u_R`

**Exact solution:**
```
u(x) = u_L + (u_R - u_L) * x/L  +  f / (2k) * x * (L - x)
```

**Test case:**
- `k = 1`, `f = 2`, `L = 1`, `u_L = 0`, `u_R = 0`
- Exact: `u(x) = x(1 - x)`
- Expected L2 convergence rate (P1): ≈ 2.0

**Verified:** `Validation.convergence.poisson_1d_convergence` (`/Validation/convergence/`) solves this exact problem with P1 elements on `n ∈ {4,8,16,32,64}`. Measured rates between the two finest meshes: L2 ≈ 2.00, H1 ≈ 1.00.

---

## 3. 2D Steady Conduction — Manufactured Solution

**Problem:** `-k ∇²u = f` on `[0,1]²`, `u = 0` on `∂Ω`

**Manufactured solution:** `u(x,y) = sin(πx) sin(πy)`

**Required source term:** `f = 2π²k sin(πx) sin(πy)`

**L2 convergence rates:**
- P1: ≈ 2.0
- P2: ≈ 3.0
- P3: ≈ 4.0

**Verified:** `Validation.convergence.poisson_2d_convergence` (h-refinement, P1) and `Validation.convergence.p_refinement_study` (P1/P2/P3 on `n ∈ {4,8,16,32}`) at `/Validation/convergence/`. Measured L2 rates between the two finest meshes: P1 ≈ 1.99, P2 ≈ 3.00, P3 ≈ 4.02.

---

## 4. 1D Convection-Diffusion — Analytical Benchmark

**Problem:** `-u'' + Pe · u' = 0` on `[0,1]`, `u(0)=0`, `u(1)=1`

**Exact solution:**
```
u(x) = (exp(Pe·x) - 1) / (exp(Pe) - 1)
```

**Test cases:**

| Pe | Expected behavior |
|----|------------------|
| 0.1 | Near-linear profile, Galerkin works fine |
| 1.0 | Slight boundary layer at x=1 |
| 10 | Thin boundary layer, Galerkin starts oscillating on coarse mesh |
| 100 | Very thin boundary layer, Galerkin completely wrong without SUPG |

**Verified:** `Validation.convergence.convection_diffusion_convergence` (`/Validation/convergence/`) solves this benchmark with both `'galerkin'` and `'supg'` for `Pe ∈ {0.1, 100}` on `n ∈ {8,16,32,64}`, reporting L2 error and the min/max nodal value. At `Pe=100`, `n=8` (`Pe_h=6.25 > 1`): standard Galerkin undershoots to `u_min ≈ -0.87` (violating the exact solution's `[0,1]` bound), while SUPG stays within `[0,1]` at every mesh size — the discrete-maximum-principle property that motivates SUPG, independent of raw L2 error.

---

## 5. Transient Heat Equation — Analytical Benchmark

**Problem:** `∂u/∂t = k ∂²u/∂x²` on `[0,1] × [0,T]`, homogeneous Dirichlet, `u₀ = sin(πx)`

**Exact solution:**
```
u(x,t) = sin(πx) · exp(-π² · k · t)
```

**Test case:** `k = 1`, `T = 0.5`, initial condition `u₀ = sin(πx)`

**Temporal convergence rates:**
- Backward Euler: ≈ 1.0
- Crank-Nicolson: ≈ 2.0

**Verified:** `TransientThermal.solver.scheme_comparison_1d` (`/TransientThermal/comparison/`) solves this benchmark on a `n=200` mesh for `dt ∈ {0.02, 0.01, 0.005}`. Measured rates: Backward Euler ≈ 0.95–0.97, Crank-Nicolson ≈ 1.80–1.95 (see `docs/science/governing_equations.md` §3 for the mesh-resolution caveat).

---

## Additional References

- **FEniCS documentation:** https://fenicsproject.org/documentation/
- **FEniCS Tutorial (Langtangen & Logg, 2016):** https://fenicsproject.org/pub/tutorial/html/
- **SUPG original paper:** Brooks, A.N. & Hughes, T.J.R. (1982). Streamline upwind/Petrov-Galerkin formulations for convection dominated flows with particular emphasis on the incompressible Navier-Stokes equations. *Comput. Methods Appl. Mech. Eng.*, 32, 199–259.
- **Ern & Guermond (2004):** Theory and Practice of Finite Elements. Springer Applied Mathematical Sciences, vol. 159.
