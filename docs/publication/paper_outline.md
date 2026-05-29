# JOSS paper.md — Pre-filled Outline

This is the working draft for the JOSS submission paper. Fill in the `[PLACEHOLDER]` sections as the project progresses. Final word count must be 750–1750 words.

Copy this file to `paper.md` at the repository root when ready to submit.

---

```markdown
---
title: 'FenicsWeb: A Web-Based Finite Element Framework for
        Convection-Diffusion Heat Transfer with Inverse Conductivity Recovery'
tags:
  - Python
  - finite element method
  - heat transfer
  - convection-diffusion
  - FEniCS
  - Django
  - inverse problems
  - SUPG stabilization
authors:
  - name: [YOUR FULL NAME]
    orcid: [YOUR ORCID — register free at orcid.org]
    affiliation: 1
affiliations:
  - name: [Your Department, Your University, Country]
    index: 1
date: [YYYY-MM-DD — date of submission]
bibliography: paper.bib
---

# Summary

FenicsWeb is an open-source, browser-based finite element solver for
heat conduction and convection-diffusion problems in one, two, and
three spatial dimensions. Built on the FEniCS finite element library
[@fenics_2015] and the Django web framework, it provides a graphical
interface for mesh generation, solver configuration, interactive
solution visualization, and solution export in VTK and XDMF formats
compatible with ParaView. The platform includes steady-state and
transient heat solvers, a convection-diffusion solver with
Streamline Upwind Petrov-Galerkin (SUPG) stabilization [@brooks_hughes_1982],
a convergence study dashboard, and an adjoint-based inverse thermal
conductivity solver powered by dolfin-adjoint [@dolfin_adjoint_2013].

# Statement of Need

Finite element simulation tools for heat transfer are widely used
in engineering and computational science, but they typically require
desktop installation of specialist software and significant programming
expertise. Existing web-based simulation platforms such as AirTrafficSim
[@airtrafficsim_2023] and OSSCAR [@osscar_2023] demonstrate demand for
browser-accessible computational tools, but none combines a research-grade
FEM backend with convection-diffusion physics, SUPG stabilization, and
inverse problem capability in a single deployable platform.

FenicsWeb fills this gap by providing:

1. A validated steady-state and transient heat solver accessible from
   any web browser, with no local installation required beyond Docker.
2. A convection-diffusion solver with SUPG stabilization and an
   interactive Peclet number study that is suitable for teaching and
   research in numerical methods.
3. An adjoint-based inverse conductivity solver that recovers the
   unknown thermal conductivity field from sparse, noisy boundary
   measurements — a capability not available in any prior web-based FEM tool.
4. A convergence study dashboard providing h-refinement and p-refinement
   error plots for independent verification of solver correctness.

The target users are researchers and students in computational mechanics,
heat transfer, and numerical methods who need rapid access to a validated
FEM solver without setting up a local simulation environment.

# Features and Functionality

**Mesh generation.** FenicsWeb supports structured 1D interval meshes,
2D rectangle meshes, and 3D box meshes via the FEniCS mesh API. Custom
meshes in Gmsh `.msh` format can be uploaded and converted automatically.

**Steady-state heat conduction.** Solves $-\nabla \cdot (k \nabla u) = f$
with Dirichlet, Neumann, and Robin boundary conditions. P1, P2, and P3
Lagrange elements are supported.

**Convection-diffusion solver.** Solves $-k\nabla^2 u + \mathbf{b} \cdot
\nabla u = f$ with standard Galerkin and SUPG-stabilized formulations.
A dedicated Peclet number study page compares both methods at Pe = 1,
10, and 100, illustrating the stabilization effect.

**Transient heat equation.** Solves $\rho c_p \partial u / \partial t =
\nabla \cdot (k \nabla u) + f$ using Backward Euler and Crank-Nicolson
time integration. Solutions are stored per timestep, viewable via a
browser timeline slider, and exportable as animated GIFs.

**Inverse conductivity solver.** Recovers the unknown conductivity
$k(\mathbf{x})$ from sparse, noisy temperature measurements by solving
a PDE-constrained optimization problem with Tikhonov regularization,
using dolfin-adjoint for automatic adjoint computation.

**Async job queue.** Long-running simulations are dispatched to a
Celery [@celery] worker queue backed by Redis, preventing HTTP request
timeouts and enabling concurrent multi-user usage.

**Solution export.** Results are exportable in VTK `.pvd`/`.vtu` and
XDMF formats for post-processing in ParaView, and as CSV for downstream
analysis.

# Verification and Validation

FenicsWeb includes a dedicated validation dashboard accessible from the
web interface. For the 1D steady-state problem with manufactured solution
$u(x) = x(1-x)$, h-refinement studies confirm L2 convergence rates of
$\mathcal{O}(h^2)$, $\mathcal{O}(h^3)$, and $\mathcal{O}(h^4)$ for
P1, P2, and P3 elements respectively, consistent with classical FEM
error estimates [@ern_guermond_2004]. For the convection-diffusion problem,
SUPG stabilization eliminates the spurious oscillations present in the
standard Galerkin solution at Pe = 100 while maintaining second-order
L2 convergence. The steady-state thermal solver was compared against the
De Vahl Davis natural convection benchmark [@de_vahl_davis_1983]; the
Nusselt number at the heated wall matches published reference values to
within [X]%.

# Research Application: Inverse Conductivity Recovery

To demonstrate research-level use, FenicsWeb was applied to the inverse
heat conduction problem in a convection-dominated domain. The study
investigated how convection strength (measured by Peclet number Pe)
affects the identifiability of the conductivity field $k(x)$ from sparse
boundary measurements.

Using a synthetic experiment with $k_{\text{true}}(x) = 1 + 0.5\sin(2\pi x)$
and 10 measurement points with Gaussian noise $\sigma = 0.01$, the
inverse solver was run at Pe = 0, 1, 5, 10, and 50. Recovery quality,
measured by L2 error $\|k - k_{\text{true}}\|_{L^2}$, degraded
monotonically with Pe: [FILL IN RESULT VALUES FROM EXPERIMENT 1].
At Pe = 50, the optimal Tikhonov regularization parameter $\alpha$
required for stable recovery was [X] times larger than at Pe = 0,
indicating that strong convection increases the ill-posedness of the
inverse problem.

A measurement sparsity study confirmed that a minimum of [N] measurement
points are required for stable recovery at Pe = 10, consistent with
the information-theoretic expectation that downstream measurements
contribute less information about upstream conductivity when Pe $\gg 1$.

These results, which are accessible interactively through the FenicsWeb
interface, represent an original numerical study not previously reported
in the web-based FEM literature.

# Acknowledgements

[PLACEHOLDER — thank supervisor, funding body, institution]

# References
```

---

## Word Count Guidance

Target sections and approximate word counts to stay within 750–1750:

| Section | Target words |
| ------- | ------------ |
| Summary | 80–100 |
| Statement of Need | 200–250 |
| Features and Functionality | 200–250 |
| Verification and Validation | 150–200 |
| Research Application | 200–250 |
| Acknowledgements | 30–50 |
| **Total** | **860–1100** |

This keeps the paper comfortably within the 1750-word ceiling while meeting the 750-word floor.

---

## Notes on Filling the Placeholders

- `[FILL IN RESULT VALUES FROM EXPERIMENT 1]` — replace after Phase 6 Experiment 1 is run; insert the actual L2 error numbers
- `[X] times larger` — replace with the actual ratio of optimal α at Pe=50 vs Pe=0
- `[N] measurement points` — replace with the minimum from Experiment 2
- ORCID — register at orcid.org, it is free and takes 2 minutes
- Acknowledgements — ask your supervisor what funding body to acknowledge
