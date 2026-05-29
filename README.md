# FenicsWeb

![CI](https://github.com/harishragulkarthik/FenicsWeb/actions/workflows/ci.yml/badge.svg)
![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)
![Python 3.10](https://img.shields.io/badge/python-3.10-blue.svg)

A web-based finite element solver for heat conduction and convection-diffusion problems, built on [FEniCS](https://fenicsproject.org) and [Django](https://djangoproject.com). Solve 1D, 2D, and 3D thermal problems directly from your browser — no local FEniCS installation required.

---

## Features

| Capability | Status |
| ---------- | ------ |
| Steady-state heat conduction (1D / 2D / 3D) | Available |
| Dirichlet, Neumann, Robin boundary conditions | In progress (Phase 1) |
| Convection-diffusion with SUPG stabilization | Planned (Phase 2) |
| Transient heat equation (Crank-Nicolson) | Planned (Phase 3) |
| h/p-refinement convergence study dashboard | Planned (Phase 4) |
| Async job queue (Celery + Redis) | Planned (Phase 5) |
| VTK / XDMF export for ParaView | Planned (Phase 5) |
| Inverse conductivity recovery (dolfin-adjoint) | Planned (Phase 6) |

---

## Installation

### Option 1 — Conda (recommended)

```bash
git clone https://github.com/harishragulkarthik/FenicsWeb.git
cd FenicsWeb
conda env create -f environment.yml
conda activate fenicsproject
python manage.py migrate
python manage.py runserver
```

Open <http://127.0.0.1:8000>

### Option 2 — Docker

```bash
git clone https://github.com/harishragulkarthik/FenicsWeb.git
cd FenicsWeb
docker-compose up --build
```

Open <http://127.0.0.1:8000>

---

## Running Tests

```bash
pytest tests/ -v
```

---

## Scientific Background

FenicsWeb uses the [Finite Element Method (FEM)](https://en.wikipedia.org/wiki/Finite_element_method) to solve partial differential equations governing heat transfer. The core equation is the steady-state Poisson equation:

```text
-∇·(k∇u) = f
```

where `u` is temperature, `k` is thermal conductivity, and `f` is a volumetric heat source. The weak formulation is discretized over a user-defined mesh using FEniCS's automated assembly pipeline.

Full mathematical details and governing equations are in [docs/science/governing_equations.md](docs/science/governing_equations.md).

---

## Project Roadmap

The full enhancement roadmap targeting PhD-application level and JOSS publication is documented in [docs/PLAN.md](docs/PLAN.md).

---

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md).

---

## License

MIT — see [LICENSE](LICENSE).

---

## Citation

If you use FenicsWeb in your research, please cite it as:

```bibtex
@software{fenics_web,
  author  = {Harish Ragul Karthik},
  title   = {FenicsWeb: A Web-Based Finite Element Framework for
             Convection-Diffusion Heat Transfer},
  year    = {2026},
  url     = {https://github.com/harishragulkarthik/FenicsWeb}
}
```

A Zenodo DOI will be added in Phase 7.
