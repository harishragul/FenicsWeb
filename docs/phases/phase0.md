# Phase 0 — Repository & CI Setup

**Duration:** 3 days
**Status:** Not Started
**Depends on:** Nothing
**Unlocks:** Everything — this is the foundation for JOSS publication eligibility

---

## Goal

Make the repository public and establish automated testing and CI from day one. This phase has no scientific content — it is entirely about starting the **6-month JOSS clock** and putting the infrastructure in place that every later phase depends on.

> JOSS explicitly rejects submissions with "a burst of commits right before submission." The 6-month public history must be genuine. Do this on Day 1.

---

## Tasks

### 0.1 Make the Repository Public

Go to GitHub → Repository Settings → Danger Zone → Change visibility → Public.

That is all. The clock starts the moment this is done.

Verify the repo URL is accessible without login before proceeding.

---

### 0.2 Add a Proper `.gitignore`

The current `.gitignore` may be missing entries. Confirm these are excluded:

```
# Python
__pycache__/
*.py[cod]
*.egg-info/
dist/
build/

# Django
*.sqlite3
/static/
/media/
/staticfiles/

# FEniCS output files
*.pvd
*.vtu
*.xdmf
*.h5
*.xml

# Environment
.env
*.env
/venv/
/env/

# IDE
.idea/
.vscode/
*.iml

# macOS
.DS_Store
```

---

### 0.3 Write the pytest Test Suite (Skeleton)

Create `tests/` at the project root. JOSS reviewers will check that tests exist and that CI runs them.

**File structure:**

```
tests/
├── __init__.py
├── conftest.py           # shared fixtures (test meshes, standard BCs)
├── test_mesh.py          # mesh generation tests
├── test_conduction.py    # steady-state conduction tests
├── test_convection.py    # placeholder — filled in Phase 2
└── test_transient.py     # placeholder — filled in Phase 3
```

**`tests/conftest.py`** — shared fixtures:
```python
import pytest
from fenics import UnitIntervalMesh, UnitSquareMesh

@pytest.fixture
def mesh_1d():
    return UnitIntervalMesh(10)

@pytest.fixture
def mesh_2d():
    return UnitSquareMesh(8, 8)
```

**`tests/test_mesh.py`** — basic mesh tests:
```python
from fenics import UnitIntervalMesh, UnitSquareMesh, UnitCubeMesh
from mesh.mesh import generate_mesh_plot

def test_interval_mesh_creates():
    mesh = UnitIntervalMesh(10)
    assert mesh.num_cells() == 10

def test_square_mesh_creates():
    mesh = UnitSquareMesh(4, 4)
    assert mesh.num_cells() > 0

def test_cube_mesh_creates():
    mesh = UnitCubeMesh(2, 2, 2)
    assert mesh.geometry().dim() == 3
```

**`tests/test_conduction.py`** — solver correctness tests:
```python
import pytest
from math import sqrt
from fenics import *
from SteadyStateThermal.conduction import solve_conduction

def test_1d_zero_source_linear_solution():
    """With f=0, Dirichlet u(0)=0, u(1)=1, solution must be linear: u(x)=x."""
    mesh = UnitIntervalMesh(20)
    u_sol, _ = solve_conduction('1D', mesh, left_bc=0, right_bc=1,
                                 top_bc=0, bottom_bc=0, front_bc=0, back_bc=0,
                                 f=0, k=1)
    # Check midpoint value
    midpoint_val = u_sol(0.5)
    assert abs(midpoint_val - 0.5) < 1e-10

def test_1d_l2_error_converges():
    """L2 error for P1 on uniform mesh must converge at O(h^2)."""
    errors = []
    for n in [8, 16, 32]:
        mesh = UnitIntervalMesh(n)
        u_sol, _ = solve_conduction('1D', mesh, left_bc=0, right_bc=0,
                                     top_bc=0, bottom_bc=0, front_bc=0, back_bc=0,
                                     f=2, k=1)
        V = u_sol.function_space()
        u_exact = interpolate(Expression("x[0]*(1-x[0])", degree=2), V)
        err = sqrt(assemble((u_sol - u_exact)**2 * dx))
        errors.append(err)
    rate = (errors[0] / errors[1])  # should be ~4 (h halves, error quarters)
    assert rate > 3.5, f"Convergence rate too low: {rate}"
```

**Placeholder files** for phases not yet implemented:
```python
# tests/test_convection.py
# Tests added in Phase 2

# tests/test_transient.py
# Tests added in Phase 3
```

---

### 0.4 Add GitHub Actions CI

**File:** `.github/workflows/ci.yml`

```yaml
name: CI

on:
  push:
    branches: [main]
  pull_request:
    branches: [main]

jobs:
  test:
    runs-on: ubuntu-latest

    steps:
      - uses: actions/checkout@v4

      - name: Set up Conda
        uses: conda-incubator/setup-miniconda@v3
        with:
          environment-file: environment.yml
          activate-environment: fenicsproject
          auto-activate-base: false

      - name: Run tests
        shell: bash -el {0}
        run: |
          conda activate fenicsproject
          pytest tests/ -v --tb=short
```

Update `environment.yml` to add `pytest`:
```yaml
name: fenicsproject
channels:
  - conda-forge
  - defaults
dependencies:
  - python=3.10
  - fenics
  - django
  - matplotlib
  - plotly
  - pytest          # ADD THIS
  - pytest-django   # ADD THIS
```

Add `pytest.ini` at the project root:
```ini
[pytest]
DJANGO_SETTINGS_MODULE = FenicsWeb.settings
python_files = tests/test_*.py
```

---

### 0.5 Add a Minimal `CONTRIBUTING.md`

JOSS requires this file. Even a short version satisfies the requirement:

```markdown
# Contributing to FenicsWeb

## Setting up the development environment
conda create -n fenicsproject -c conda-forge fenics django matplotlib plotly pytest
conda activate fenicsproject

## Running tests
pytest tests/ -v

## Reporting issues
Open an issue on GitHub with a description of the problem and steps to reproduce.

## Submitting changes
1. Fork the repository
2. Create a branch: git checkout -b feature/your-feature
3. Make changes and add tests
4. Ensure tests pass: pytest tests/
5. Open a pull request
```

---

### 0.6 Overhaul the README

The current README is 4 lines. JOSS reviewers and PhD guides look at the README first.

Minimum sections required:
- Project description (1–2 sentences)
- Screenshot of the running application
- Installation instructions (conda + pip paths)
- Quickstart (how to run, URL to open)
- Feature list
- Scientific background (1 paragraph on FEM and what problems are solved)
- License + citation badge

The full README overhaul is in Phase 5. Phase 0 just needs enough to not embarrass the project when the repo goes public.

---

## Files to Create / Modify

| File | Action |
| ---- | ------ |
| `.gitignore` | Update with FEniCS + Django entries |
| `tests/__init__.py` | New — empty |
| `tests/conftest.py` | New — shared fixtures |
| `tests/test_mesh.py` | New — mesh tests |
| `tests/test_conduction.py` | New — solver correctness tests |
| `tests/test_convection.py` | New — placeholder |
| `tests/test_transient.py` | New — placeholder |
| `.github/workflows/ci.yml` | New — GitHub Actions |
| `environment.yml` | Add pytest, pytest-django |
| `pytest.ini` | New — pytest configuration |
| `CONTRIBUTING.md` | New — minimal contributor guide |

---

## Definition of Done

- [ ] Repository is public on GitHub
- [ ] `pytest tests/` runs without errors locally
- [ ] GitHub Actions CI badge shows green on the README
- [ ] `CONTRIBUTING.md` exists
- [ ] `.gitignore` covers all FEniCS/Django/Python artifacts
- [ ] At least 3 passing tests in the suite
