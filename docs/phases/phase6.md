# Phase 6 — Documentation & Testing

**Duration:** 2 weeks
**Status:** Not Started
**Depends on:** Phases 1–5 (all solvers, verification framework, and export/import implemented)
**Unlocks:** Phase 7 (publication readiness)

---

## Goal

Implement Goals 10 and 11: a complete set of scientific documentation
(governing equations, weak forms, numerical methods, verification,
benchmarks, examples, theory notes) and an automated test suite reaching a
minimum of 90% coverage across mesh generation, boundary conditions, the
steady solver, the transient solver, the convection-diffusion solver, and
the export modules.

---

## Tasks

### 6.1 Scientific Documentation (Goal 10)

`docs/science/governing_equations.md` and `docs/science/benchmark_references.md`
already exist from Phase 0/1. Expand `docs/science/` to cover all required
sections:

| File | Required content |
| ---- | ----------------- |
| `governing_equations.md` | Steady conduction, convection-diffusion, transient heat equation — strong forms, boundary conditions, units |
| `weak_forms.md` | Full weak-form derivations: Dirichlet/Neumann/Robin (Phase 1), Galerkin and SUPG convection-diffusion (Phase 2), Backward Euler and Crank-Nicolson time discretization (Phase 3) |
| `numerical_methods.md` | FE spaces (P1/P2/P3), time-integration schemes, SUPG stabilization parameter `τ` and its derivation |
| `verification.md` | L2/H1 error norm definitions, convergence-rate estimation methodology, how to reproduce the Phase 4 convergence tables/plots |
| `benchmark_references.md` | 1D/2D conduction exact solutions, De Vahl Davis benchmark — expand with the Phase 4 results table |
| `theory_notes.md` | Peclet number and its physical meaning, numerical diffusion, stability/conditioning notes for convection-dominated problems |

Add `docs/examples/` with one short workflow document per solver (steady
conduction, convection-diffusion, transient, verification study, mesh
import + export), each showing the inputs used and the expected output
(plot/table), so results are reproducible from the documentation alone.

---

### 6.2 README Overhaul

Replace the current minimal README with sections covering:

- **Badges row** — CI status, license, Python version (JOSS badge added in Phase 7)
- **Screenshot** — actual browser screenshot of a solver result
- **What it does** — 1–2 sentences on the tool and its physics scope
- **Installation** — conda/pip instructions with copy-pasteable commands
- **Quickstart** — how to run and open the first solver in a few steps
- **Features table** — solvers, BC types, FE orders, export formats, mesh import
- **Scientific background** — short paragraph on the FEM and equations solved, linking to `docs/science/governing_equations.md`
- **Verification** — short paragraph + link to `docs/science/verification.md` and the Phase 4 convergence results
- **Citation** — BibTeX block placeholder, completed with the Zenodo DOI in Phase 7
- **License** — MIT

---

### 6.3 API Docstrings

Add a one-line docstring to every public function across the solver
modules. Standard format:

```python
def solve_conduction(dimension, mesh, left_bc, right_bc, ...):
    """Solve steady-state heat conduction via FEniCS FEM; return (u_sol, solution_list)."""
```

Functions to cover:

- `SteadyStateThermal/conduction.py` — `solve_conduction`, `generate_solution_plot`
- `ConvectionDiffusion/solver.py` — Galerkin and SUPG solve functions, Peclet number helper
- `TransientThermal/solver.py` — Backward Euler / Crank-Nicolson solve functions
- `mesh/mesh.py`, `mesh/export.py`, `mesh/gmsh_import.py`
- `Validation/convergence.py`, `Validation/benchmarks.py`

---

### 6.4 `CONTRIBUTING.md` (full upgrade from Phase 0 skeleton)

Expand to include:

- Development environment setup
- How to run tests (`pytest tests/ -v --cov`)
- How to add a new solver, pointing to existing solvers as templates
- Code style guide (PEP 8, max line length 100)
- How to submit a pull request

---

### 6.5 Test Coverage to 90% (Goal 11)

- Add `pytest-cov` and a coverage configuration (`.coveragerc` or
  `pyproject.toml` `[tool.coverage]` section) with a 90% minimum threshold.
- Replace the placeholder `tests/test_convection.py` and
  `tests/test_transient.py` with real tests for the Phase 2 and Phase 3
  solvers (correctness vs. analytical solutions, SUPG vs. Galerkin
  comparison, Backward Euler / Crank-Nicolson convergence).
- Add `tests/test_verification.py` for the Phase 4 convergence-rate
  estimation code, and `tests/test_export.py` / `tests/test_mesh_import.py`
  for the Phase 5 export/import modules (these may already exist from
  Phase 5; ensure they contribute to the coverage target).
- Add a CI step that runs `pytest --cov --cov-fail-under=90` so coverage
  regressions break the build.

---

## Files to Create

| File | Purpose |
| ---- | ------- |
| `docs/science/weak_forms.md` | Full weak-form derivations for all BCs, convection-diffusion, transient schemes |
| `docs/science/numerical_methods.md` | FE spaces, time integration, SUPG stabilization |
| `docs/science/verification.md` | Error norms + convergence methodology |
| `docs/science/theory_notes.md` | Peclet number, numerical diffusion, stability |
| `docs/examples/*.md` | One reproducible example workflow per solver |
| `docs/usage/getting_started.md` | First-solve walkthrough |
| `docs/usage/solvers.md` | Solver reference (inputs/outputs per solver) |
| `.coveragerc` | Coverage configuration with 90% threshold |
| `CONTRIBUTING.md` | Full contributor guide |

**Files to modify:**

| File | Change |
| ---- | ------ |
| `README.md` | Full overhaul per 6.2 |
| `.github/workflows/ci.yml` | Add `pytest --cov --cov-fail-under=90` step |
| `tests/test_convection.py`, `tests/test_transient.py` | Replace placeholders with real tests |
| All solver Python files | Add one-line docstrings to public functions |

---

## Definition of Done

- [ ] `docs/science/` covers governing equations, weak forms, numerical methods, verification, benchmarks, examples, and theory notes
- [ ] README has badges, screenshot, installation, quickstart, features table, scientific background, verification section, citation placeholder, license
- [ ] All public solver functions have one-line docstrings
- [ ] `CONTRIBUTING.md` covers setup, tests, adding a solver, PR process
- [ ] `docs/usage/` has getting-started and solver-reference pages
- [ ] `pytest --cov` reports ≥ 90% coverage across mesh, BCs, steady/transient/convection solvers, and export modules
- [ ] CI fails if coverage drops below 90%
