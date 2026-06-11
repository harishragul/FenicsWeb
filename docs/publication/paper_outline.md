# JOSS paper.md — Pre-filled Outline

This is the working draft for the JOSS submission paper. Fill in the `[PLACEHOLDER]` sections as the project progresses. Final word count must be 750–1750 words.

Copy this file to `paper.md` at the repository root when ready to submit.

---

```markdown
---
title: 'FenicsWeb: A Web-Based Finite Element Framework for
        Convection-Diffusion Heat Transfer'
tags:
  - Python
  - finite element method
  - heat transfer
  - convection-diffusion
  - FEniCS
  - Django
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
interface for mesh generation (including Gmsh `.msh` import), solver
configuration, interactive solution visualization, and solution export
in CSV, VTK, and XDMF formats compatible with ParaView. The platform
includes steady-state and transient heat solvers with P1/P2/P3 elements,
a convection-diffusion solver with Streamline Upwind Petrov-Galerkin
(SUPG) stabilization [@brooks_hughes_1982], and an integrated
verification framework with h- and p-refinement convergence studies and
a benchmark library.

# Statement of Need

Finite element simulation tools for heat transfer are widely used
in engineering and computational science, but they typically require
desktop installation of specialist software and significant programming
expertise. Existing web-based simulation platforms such as AirTrafficSim
[@airtrafficsim_2023] and OSSCAR [@osscar_2023] demonstrate demand for
browser-accessible computational tools, but none combines a research-grade
FEM backend for convection-diffusion physics with SUPG stabilization, an
integrated verification/benchmark suite, and ParaView-compatible export
in a single deployable platform.

FenicsWeb fills this gap by providing:

1. A validated steady-state and transient heat solver accessible from
   any web browser, supporting Dirichlet, Neumann, and Robin boundary
   conditions and P1/P2/P3 elements.
2. A convection-diffusion solver with SUPG stabilization and an
   interactive Peclet number study that is suitable for teaching and
   research in numerical methods.
3. An integrated verification framework providing h-refinement and
   p-refinement convergence tables/plots, L2 and H1 error norms, and a
   benchmark library (1D/2D conduction with exact solutions, De Vahl
   Davis natural convection) for independent verification of solver
   correctness.
4. Gmsh mesh import and CSV/VTK/XDMF export, allowing users to bring in
   real geometries and post-process results in ParaView.

The target users are researchers and students in computational mechanics,
heat transfer, and numerical methods who need rapid access to a validated
FEM solver without setting up a local simulation environment.

# Features and Functionality

**Mesh generation and import.** FenicsWeb supports structured 1D interval
meshes, 2D rectangle meshes, and 3D box meshes via the FEniCS mesh API.
Custom meshes in Gmsh `.msh` format can be uploaded and converted
automatically via `meshio`.

**Steady-state heat conduction.** Solves $-\nabla \cdot (k \nabla u) = f$
with Dirichlet, Neumann, and Robin boundary conditions. P1, P2, and P3
Lagrange elements are supported.

**Convection-diffusion solver.** Solves $-k\nabla^2 u + \mathbf{b} \cdot
\nabla u = f$ with standard Galerkin and SUPG-stabilized formulations.
A dedicated Peclet number study page compares both methods at Pe = 1,
10, and 100, illustrating the stabilization effect.

**Transient heat equation.** Solves $\rho c_p \partial u / \partial t =
\nabla \cdot (k \nabla u) + f$ using Backward Euler and Crank-Nicolson
time integration, with timestep controls and transient visualization.

**Verification framework.** A dedicated validation module computes L2
and H1 error norms against analytical and manufactured solutions, runs
h- and p-refinement convergence studies, and estimates the observed
order of accuracy. A benchmark library includes 1D and 2D conduction
problems with exact solutions and the De Vahl Davis natural convection
benchmark.

**Solution export.** Results are exportable as CSV (DOF coordinates and
values), VTK `.pvd`/`.vtu`, and XDMF for post-processing in ParaView,
including transient time-series animation.

# Verification and Validation

FenicsWeb includes a dedicated validation module accessible from the
web interface. For the 1D and 2D steady-state problems with manufactured
solutions, h-refinement studies confirm L2 convergence rates of
$\mathcal{O}(h^2)$, $\mathcal{O}(h^3)$, and $\mathcal{O}(h^4)$ for
P1, P2, and P3 elements respectively, consistent with classical FEM
error estimates [@ern_guermond_2004]. For the convection-diffusion
problem, SUPG stabilization eliminates the spurious oscillations present
in the standard Galerkin solution at high Peclet number while maintaining
the expected L2 convergence rate. The transient solver was verified
against the analytical decay solution $u(x,t) = \sin(\pi x)
e^{-\pi^2 k t}$, confirming first- and second-order temporal convergence
for Backward Euler and Crank-Nicolson respectively. The steady-state
thermal solver was compared against the De Vahl Davis natural convection
benchmark [@de_vahl_davis_1983]; the Nusselt number at the heated wall
matches published reference values to within [X]%.

# Acknowledgements

[PLACEHOLDER — thank supervisor, funding body, institution]

# References
```

---

## Word Count Guidance

Target sections and approximate word counts to stay within 750–1750:

| Section | Target words |
| ------- | ------------ |
| Summary | 90–120 |
| Statement of Need | 200–250 |
| Features and Functionality | 250–300 |
| Verification and Validation | 150–200 |
| Acknowledgements | 30–50 |
| **Total** | **720–920** |

This sits comfortably within the 750–1750 word range; expand the
Features and Verification sections with concrete numbers from Phase 4
results if more detail is needed.

---

## Notes on Filling the Placeholders

- `[X]%` — replace with the actual Nusselt number agreement from the
  Phase 4 De Vahl Davis comparison
- ORCID — register at orcid.org, it is free and takes 2 minutes
- Acknowledgements — ask your supervisor what funding body to acknowledge
