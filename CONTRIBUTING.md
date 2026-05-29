# Contributing to FenicsWeb

Thank you for your interest in contributing to FenicsWeb.

## Setting Up the Development Environment

### Option 1 — Conda (recommended)

```bash
conda env create -f environment.yml
conda activate fenicsproject
```

### Option 2 — Docker

```bash
docker-compose up --build
```

## Running the Tests

```bash
pytest tests/ -v
```

All tests must pass before submitting a pull request.

## Running the Development Server

```bash
python manage.py migrate
python manage.py runserver
```

Open http://127.0.0.1:8000 in your browser.

## Project Structure

```
FenicsWeb/
├── mesh/                  # Mesh generation
├── SteadyStateThermal/    # Steady-state conduction solver
├── ConvectionDiffusion/   # Convection-diffusion solver (Phase 2)
├── TransientThermal/      # Transient heat solver (Phase 3)
├── Validation/            # Convergence studies and benchmarks (Phase 4)
├── Research/              # Inverse problem solver (Phase 6)
├── tests/                 # pytest test suite
└── docs/                  # Project documentation and plan
```

## Adding a New Solver

1. Create a new Django app: `python manage.py startapp MySolver`
2. Add a `solver.py` with the FEniCS solve function
3. Add a `forms.py` for the user input form
4. Add a `views.py`, `urls.py`, and a template under `templates/`
5. Register the app in `FenicsWeb/settings.py` and `FenicsWeb/urls.py`
6. Add tests to `tests/test_mysolverapp.py` — at minimum one correctness test and one convergence test

## Code Style

- PEP 8, maximum line length 100
- One-line docstrings on all public functions
- No commented-out code in pull requests

## Reporting Issues

Open an issue on GitHub with:
- A short description of the problem
- Steps to reproduce
- Expected vs actual behaviour
- Python/FEniCS version (`python --version`, `python -c "import fenics; print(fenics.__version__)"`)

## Submitting Changes

1. Fork the repository
2. Create a branch: `git checkout -b feature/your-feature`
3. Make changes and add tests
4. Confirm tests pass: `pytest tests/ -v`
5. Open a pull request with a clear description of what changed and why
