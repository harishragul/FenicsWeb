# FenicsWeb — Scientific Computing Enhancement Plan

**Project:** Web-based Finite Element solver for heat conduction and convection
**Stack:** FEniCS · Django · SQLite · Plotly · Matplotlib
**Goal:** Elevate FenicsWeb from a student demo to a compact, scientifically
rigorous finite element research platform suitable for a master's thesis
enhancement, a PhD application portfolio, and a JOSS-style software
publication.

**Publication Target (primary):**

> Journal of Open Source Software (JOSS) — peer-reviewed, DOI-bearing, Scopus/WoS indexed

---

## CRITICAL: Start the 6-Month Clock Today

JOSS requires a public repository with **at least 6 months of genuine commit history** before submission. The clock starts the day the repo is made public. Make the repo public on GitHub immediately — even if the code is not perfect yet. Commit regularly throughout all phases.

> **Action:** Go to GitHub repo settings → change visibility to Public → done.

---

## Status Tracker

| Phase | Name | Status | Goals Covered | Publication Impact |
| ----- | ---- | ------ | -------------- | ------------------- |
| [Phase 0](phases/phase0.md) | Repo & CI Setup | **Complete** | — | **Starts 6-month JOSS clock** |
| [Phase 1](phases/phase1.md) | Foundation Repair | **Complete** | Goal 1 (BCs) | First validated result |
| [Phase 2](phases/phase2.md) | Convection-Diffusion + SUPG | Not Started | Goals 2, 3 | Core scientific contribution |
| [Phase 3](phases/phase3.md) | Transient Heat Equation | Not Started | Goal 4 | Broadens solver scope |
| [Phase 4](phases/phase4.md) | Higher-Order Elements, Verification & Benchmarks | Not Started | Goals 5, 6, 7 | V&V section of paper |
| [Phase 5](phases/phase5.md) | Scientific Data Export & Mesh Import | Not Started | Goals 8, 9 | ParaView-compatible reproducibility |
| [Phase 6](phases/phase6.md) | Documentation & Testing | Not Started | Goals 10, 11 | JOSS software requirements |
| [Phase 7](phases/phase7.md) | Publication Preparation | Not Started | Goal 12 | JOSS submission |

**Total estimated duration:** ~21 weeks

---

## Publication Strategy

```text
Week 1      → Make repo public (Phase 0) — 6-month JOSS clock starts
Week 1–3    → Phase 1: Fix bugs, add BCs, first tests
Week 4–7    → Phase 2: Convection-diffusion + SUPG
Week 8–11   → Phase 3: Transient solver
Week 12–14  → Phase 4: Higher-order elements, verification framework, benchmarks
Week 15–17  → Phase 5: CSV/VTK/XDMF export, Gmsh mesh import
Week 18–19  → Phase 6: Documentation pass + test coverage to 90%
Week 20–21  → Phase 7: Zenodo + JOSS submission
```

By the end of Phase 7 the repository will have ~21 weeks of public commit history, approaching the 6-month JOSS requirement.

---

## Architecture Overview

```text
FenicsWeb/
├── mesh/                        # Mesh generation (interval, rectangle, box) + Gmsh import (Phase 5)
├── SteadyStateThermal/          # Steady-state conduction (Phase 1 — complete)
├── ConvectionDiffusion/         # NEW — Phase 2 (convection-diffusion + SUPG)
├── TransientThermal/            # NEW — Phase 3
├── Validation/                  # NEW — Phase 4 (verification framework, benchmark library)
├── tests/                       # Phase 0 — pytest suite, target 90% coverage (Goal 11)
│   ├── test_mesh.py
│   ├── test_conduction.py
│   ├── test_convection.py
│   ├── test_transient.py
│   ├── test_verification.py     # Phase 4
│   └── test_export.py           # Phase 5
├── .github/
│   └── workflows/
│       └── ci.yml               # Phase 0 — GitHub Actions CI
├── CONTRIBUTING.md              # Phase 0 skeleton; full upgrade in Phase 6
└── docs/
    ├── PLAN.md                  # This file
    ├── phases/
    │   ├── phase0.md
    │   ├── phase1.md
    │   ├── phase2.md
    │   ├── phase3.md
    │   ├── phase4.md
    │   ├── phase5.md
    │   ├── phase6.md
    │   └── phase7.md
    ├── science/
    │   ├── governing_equations.md
    │   ├── benchmark_references.md
    │   └── weak_forms.md
    └── publication/
        ├── joss_checklist.md
        └── paper_outline.md
```

---

## Governing Equations Summary

| Problem | PDE | Phase |
| ------- | --- | ----- |
| Steady conduction | `-∇·(k∇u) = f` | Baseline |
| Steady convection-diffusion | `-k∇²u + b·∇u = f` | Phase 2 |
| Transient heat equation | `ρcp ∂u/∂t = ∇·(k∇u) + f` | Phase 3 |

---

## Key Design Decisions

### Boundary Condition Types

- **Dirichlet** (`u = g`): prescribed temperature — implemented in Phase 1
- **Neumann** (`-k ∂u/∂n = q`): prescribed heat flux — implemented in Phase 1
- **Robin** (`-k ∂u/∂n = h(u - u∞)`): convective cooling — implemented in Phase 1

### Finite Element Spaces

- **P1** (linear Lagrange): baseline, currently used
- **P2/P3** (higher-order): added in Phase 4 for p-refinement study (Goal 5)

### Time Integration (Phase 3)

- **Backward Euler**: 1st-order, unconditionally stable
- **Crank-Nicolson**: 2nd-order, unconditionally stable — primary scheme

### Stabilization (Phase 2)

- **SUPG** (Streamline Upwind Petrov-Galerkin): standard stabilization for convection-dominated problems
- Activates when cell Peclet number `Pe_h = |b|h/(2k) > 1`

---

## Scientific Contributions Checklist

- [x] Neumann and Robin boundary conditions — Phase 1 complete
- [x] Analytical validation table (1D exact solutions) — Phase 1 complete
- [ ] Convection-diffusion solver with Peclet number computation — Phase 2 (Goal 2)
- [ ] SUPG stabilization with stabilized vs. unstabilized comparison — Phase 2 (Goal 3)
- [ ] Transient solver with two time-stepping schemes + temporal convergence — Phase 3 (Goal 4)
- [ ] P1/P2/P3 solution and convergence comparison — Phase 4 (Goal 5)
- [ ] L2/H1 error norms, mesh-refinement convergence tables/plots, order-of-accuracy estimates — Phase 4 (Goal 6)
- [ ] 1D/2D conduction benchmarks vs. exact solutions + De Vahl Davis benchmark — Phase 4 (Goal 7)

## Software Contributions Checklist

- [x] Fixed solution visualization (temperature colormap, not mesh replot) — Phase 1 complete
- [x] Interactive Plotly 3D solution field with temperature colouring — Phase 1 complete
- [ ] CSV / VTK / XDMF export, ParaView-compatible — Phase 5 (Goal 8)
- [ ] Gmsh `.msh` mesh import — Phase 5 (Goal 9)
- [ ] Scientific documentation: governing equations, weak forms, numerical methods, verification, benchmarks, examples, theory notes — Phase 6 (Goal 10)
- [ ] Automated test suite ≥ 90% coverage across mesh, BCs, steady/transient/convection solvers, export — Phase 6 (Goal 11)

## Publication Readiness Checklist

- [x] Repo public on GitHub — Phase 0 complete (6-month clock running)
- [x] GitHub Actions CI passing — Phase 0 complete
- [x] pytest suite passing — Phase 0/1 complete
- [x] CONTRIBUTING.md skeleton — Phase 0 complete
- [ ] README: full overhaul with screenshots, citation block — Phase 6
- [ ] CONTRIBUTING.md: full upgrade — Phase 6
- [ ] CITATION.cff for software citation — Phase 7
- [ ] Zenodo DOI for software archive — Phase 7
- [ ] JOSS paper.md (750–1750 words) drafted — Phase 7
- [ ] JOSS submission submitted — Phase 7

---

## Precedents (Tools Published in JOSS Using FEniCS)

| Tool | JOSS Issue | Description |
| ---- | ---------- | ----------- |
| dolfin-adjoint | joss.01292 | Automated adjoints for FEniCS |
| FEniCS-arclength | joss.05727 | Numerical continuation for solid mechanics |
| ADIOS4DOLFINx | joss.06451 | Checkpointing framework |
| fenicsx-beat | joss.08416 | Cardiac electrophysiology |
| pulse | joss.01539 | Cardiac mechanics |
| AirTrafficSim | joss.04916 | Web-based simulation platform (Flask) — closest structural precedent |

FenicsWeb is more technically serious than AirTrafficSim and follows an established JOSS publication pattern.

---

## References

- De Vahl Davis, G. (1983). *Int. J. Num. Methods Fluids*, 3, 249–264.
- Langtangen, H.P. & Logg, A. (2016). *Solving PDEs in Python: The FEniCS Tutorial I.* Springer.
- Brooks, A.N. & Hughes, T.J.R. (1982). *Comput. Methods Appl. Mech. Eng.*, 32, 199–259.
- Ern, A. & Guermond, J-L. (2004). *Theory and Practice of Finite Elements.* Springer.
