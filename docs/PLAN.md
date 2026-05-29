# FenicsWeb — PhD Application & Publication Enhancement Plan

**Project:** Web-based Finite Element solver for heat conduction and convection  
**Stack:** FEniCS · Django · SQLite · Plotly · Matplotlib · Docker  
**Goal:** Elevate from a student demo to a thesis-worthy, peer-reviewed-publishable computational research platform  

**Thesis Title (target):**
> *"A Web-Based Finite Element Framework for Convection-Diffusion Heat Transfer: Verification, Validation, and Inverse Conductivity Recovery"*

**Publication Target (primary):**
> Journal of Open Source Software (JOSS) — peer-reviewed, DOI-bearing, Scopus/WoS indexed

**Preprint Target:**
> arXiv cs.CE (Computational Engineering) + math.NA (Numerical Analysis) cross-list

---

## CRITICAL: Start the 6-Month Clock Today

JOSS requires a public repository with **at least 6 months of genuine commit history** before submission. The clock starts the day the repo is made public. Make the repo public on GitHub immediately — even if the code is not perfect yet. Commit regularly throughout all phases.

> **Action:** Go to GitHub repo settings → change visibility to Public → done.

---

## Status Tracker

| Phase | Name | Status | Duration | Publication Impact |
| ----- | ---- | ------ | -------- | ------------------ |
| [Phase 0](phases/phase0.md) | Repo & CI Setup | **Complete** | 3 days | **Starts 6-month JOSS clock** |
| [Phase 1](phases/phase1.md) | Foundation Repair | Not Started | 2 weeks | First validated result |
| [Phase 2](phases/phase2.md) | Convection-Diffusion + SUPG | Not Started | 4 weeks | Core scientific contribution |
| [Phase 3](phases/phase3.md) | Transient Heat Equation | Not Started | 4 weeks | Broadens solver scope |
| [Phase 4](phases/phase4.md) | Scientific Rigor & Benchmarks | Not Started | 3 weeks | V&V section of paper |
| [Phase 5](phases/phase5.md) | Platform Maturity + Docs | Not Started | 3 weeks | JOSS software requirements |
| [Phase 6](phases/phase6.md) | Research Contribution | Not Started | 4 weeks | Novel result for thesis/paper |
| [Phase 7](phases/phase7.md) | Publication Preparation | Not Started | 2 weeks | JOSS + arXiv submission |

**Total estimated duration:** ~23 weeks (~6 months — aligned with JOSS clock)

---

## Publication Strategy

```
Week 1      → Make repo public (Phase 0) — 6-month JOSS clock starts
Week 1–3    → Phase 1: Fix bugs, add BCs, first tests
Week 4–7    → Phase 2: Convection-diffusion + SUPG
Week 8–11   → Phase 3: Transient solver
Week 12–14  → Phase 4: Convergence + benchmarks
Week 15–17  → Phase 5: Async, docs, VTK export
Week 18–21  → Phase 6: Inverse problem research
Week 22–23  → Phase 7: arXiv preprint + JOSS submission
```

arXiv can be posted once Phase 4 is done (week 14). JOSS submission follows Phase 7 (week 23), by which time the repo will have 23 weeks of public history.

---

## Architecture Overview

```
FenicsWeb/
├── mesh/                        # Mesh generation (interval, rectangle, box)
├── SteadyStateThermal/          # Steady-state conduction (Phase 1 target)
├── ConvectionDiffusion/         # NEW — Phase 2
├── TransientThermal/            # NEW — Phase 3
├── Validation/                  # NEW — Phase 4 (convergence, benchmarks)
├── Jobs/                        # NEW — Phase 5 (Celery async jobs)
├── Research/                    # NEW — Phase 6 (inverse problem)
├── tests/                       # NEW — Phase 0 (pytest test suite)
│   ├── test_mesh.py
│   ├── test_conduction.py
│   ├── test_convection.py
│   └── test_transient.py
├── .github/
│   └── workflows/
│       └── ci.yml               # NEW — Phase 0 (GitHub Actions)
├── CONTRIBUTING.md              # NEW — Phase 5 (JOSS requirement)
└── docs/
    ├── PLAN.md                  # This file
    ├── phases/
    │   ├── phase0.md            # NEW — Repo & CI setup
    │   ├── phase1.md
    │   ├── phase2.md
    │   ├── phase3.md
    │   ├── phase4.md
    │   ├── phase5.md
    │   ├── phase6.md
    │   └── phase7.md            # NEW — Publication prep
    ├── science/
    │   ├── governing_equations.md
    │   ├── benchmark_references.md
    │   └── weak_forms.md
    └── publication/
        ├── joss_checklist.md    # NEW — JOSS pass/fail requirements
        └── paper_outline.md     # NEW — JOSS paper.md template
```

---

## Governing Equations Summary

| Problem | PDE | Phase |
|---------|-----|-------|
| Steady conduction | `-∇·(k∇u) = f` | Baseline |
| Steady convection-diffusion | `-k∇²u + b·∇u = f` | Phase 2 |
| Transient heat equation | `ρcp ∂u/∂t = ∇·(k∇u) + f` | Phase 3 |
| Inverse heat conduction (novel) | Recover `k(x)` from sparse noisy measurements under convection | Phase 6 |

---

## Key Design Decisions

### Boundary Condition Types
- **Dirichlet** (`u = g`): prescribed temperature — already implemented
- **Neumann** (`-k ∂u/∂n = q`): prescribed heat flux — Phase 1 addition
- **Robin** (`-k ∂u/∂n = h(u - u∞)`): convective cooling — Phase 1 addition; also the inverse problem domain

### Finite Element Spaces
- **P1** (linear Lagrange): baseline, currently used
- **P2/P3** (higher-order): added in Phase 4 for p-refinement study

### Time Integration (Phase 3)
- **Backward Euler**: 1st-order, unconditionally stable
- **Crank-Nicolson**: 2nd-order, unconditionally stable — primary scheme

### Stabilization (Phase 2)
- **SUPG** (Streamline Upwind Petrov-Galerkin): standard stabilization for convection-dominated problems
- Activates when cell Peclet number `Pe_h = |b|h/(2k) > 1`

---

## Scientific Contributions Checklist

- [ ] Neumann and Robin boundary conditions
- [ ] SUPG stabilization with Peclet number analysis
- [ ] h-convergence rates (L2 and H1 norms, P1/P2/P3)
- [ ] Transient solver with two time-stepping schemes + temporal convergence
- [ ] De Vahl Davis benchmark comparison
- [ ] Analytical validation table (1D exact solutions)
- [ ] **Novel:** Inverse conductivity recovery in convection-dominated domain with sparse noisy measurements — identifiability + regularization study (Phase 6)

## Software Contributions Checklist

- [ ] Fixed solution visualization (temperature colormap, not mesh replot)
- [ ] Interactive Plotly 3D solution field
- [ ] Async job queue (Celery + Redis)
- [ ] VTK / XDMF export for ParaView
- [ ] Custom mesh upload (Gmsh `.msh`)
- [ ] Parameter sweep / sensitivity study runner
- [ ] User accounts and job history
- [ ] Convergence study dashboard

## Publication Readiness Checklist

- [ ] Repo public on GitHub (Phase 0)
- [ ] pytest test suite covering all solvers (Phase 0–3)
- [ ] GitHub Actions CI passing (Phase 0)
- [ ] README: installation, usage, examples with screenshots (Phase 5)
- [ ] CONTRIBUTING.md (Phase 5)
- [ ] Zenodo DOI for software archive (Phase 7)
- [ ] arXiv preprint posted (Phase 7)
- [ ] JOSS paper.md (750–1750 words) drafted (Phase 7)
- [ ] JOSS submission submitted (Phase 7)

---

## Precedents (Tools Published in JOSS Using FEniCS)

| Tool | JOSS Issue | Description |
| ---- | --------- | ----------- |
| dolfin-adjoint | joss.01292 | Automated adjoints for FEniCS — directly used in Phase 6 |
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
- Farrell, P.E. et al. (2013). dolfin-adjoint. *SIAM J. Sci. Comput.*, 35(4), C369–C393.
- Ern, A. & Guermond, J-L. (2004). *Theory and Practice of Finite Elements.* Springer.
