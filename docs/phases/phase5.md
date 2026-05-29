# Phase 5 — Platform Maturity + Documentation

**Duration:** 3 weeks
**Status:** Not Started
**Depends on:** Phase 1, 2, 3 (all solvers working)
**Unlocks:** Phase 6 (research features need async infrastructure)

---

## Goal

Transform the project from a demo into a production-quality research platform. These additions make the tool usable by other researchers — and the documentation tasks in 5.7 are non-negotiable requirements for JOSS submission.

---

## Tasks

### 5.1 Async Job Queue (Celery + Redis)

**Problem:** A 3D solve or a transient simulation can take minutes. Synchronous Django views will time out (default: 30s). Any serious simulation must be async.

**Stack:** Celery (Python distributed task queue) + Redis (message broker).

**Installation:**

```bash
conda install -c conda-forge celery redis-py
```

**Docker addition** — add to `docker-compose.yml`:

```yaml
services:
  redis:
    image: redis:7-alpine
    ports:
      - "6379:6379"

  celery:
    build: .
    command: celery -A FenicsWeb worker --loglevel=info
    depends_on:
      - redis
    volumes:
      - .:/app
```

**`Jobs/tasks.py` skeleton:**

```python
from celery import shared_task

@shared_task
def run_steady_conduction(job_id):
    job = SolverJob.objects.get(id=job_id)
    job.status = 'running'
    job.save()
    try:
        # ... run solver ...
        job.status = 'completed'
        job.result_json = json.dumps(result)
    except Exception as e:
        job.status = 'failed'
        job.error = str(e)
    job.save()
```

**View flow:**

1. User submits form → create `SolverJob` record → enqueue Celery task → return job ID
2. Redirect to `/jobs/<job_id>/status/`
3. Status page polls `/jobs/<job_id>/api/` every 2 seconds via JavaScript
4. When status = completed, redirect to results page

---

### 5.2 Jobs Dashboard

**`Jobs/models.py`:**

```python
class SolverJob(models.Model):
    STATUS = [
        ('queued', 'Queued'),
        ('running', 'Running'),
        ('completed', 'Completed'),
        ('failed', 'Failed'),
    ]
    user = models.ForeignKey(User, on_delete=models.CASCADE)
    solver_type = models.CharField(max_length=50)
    status = models.CharField(max_length=20, choices=STATUS, default='queued')
    created_at = models.DateTimeField(auto_now_add=True)
    completed_at = models.DateTimeField(null=True)
    params_json = models.TextField()
    result_json = models.TextField(null=True)
    error = models.TextField(null=True)
```

Dashboard page: table of all user jobs with status, solver type, elapsed time, and link to results. Allows re-running or deleting old jobs.

---

### 5.3 VTK / XDMF Solution Export

VTK (`.vtu`) and XDMF (`.xdmf`) are the standard formats for ParaView — the tool researchers use for publication-quality post-processing.

**Implementation:**

```python
from fenics import XDMFFile, File

def export_vtk(u_sol, filename):
    vtkfile = File(f'static/exports/{filename}.pvd')
    vtkfile << u_sol

def export_xdmf(u_sol, mesh, filename):
    with XDMFFile(f'static/exports/{filename}.xdmf') as f:
        f.write(mesh)
        f.write(u_sol)
```

Add "Download as VTK" and "Download as XDMF" buttons to every results page.

---

### 5.4 Custom Mesh Upload (Gmsh)

The `CustomMesh` form class is already stubbed out (commented) in `mesh/forms.py:62`.

**Implementation steps:**

- Accept `.msh` (Gmsh format 2 or 4) uploads via `forms.FileField`
- Convert with `dolfin-convert` or `meshio`:

```bash
dolfin-convert mesh.msh mesh.xml
```

- Store the converted `.xml` in `media/meshes/` and reference it in the `Mesh` model

Custom meshes allow solving on real engineering geometries (heat sinks, PCB boards, building sections), which makes the tool genuinely useful beyond textbook problems.

---

### 5.5 Parameter Sweep / Sensitivity Study

Allow users to sweep any scalar parameter over a range and run all cases automatically.

**`ParameterSweepForm` fields:**

- `sweep_param` — which parameter to sweep (k, f, bx, by, dt)
- `sweep_min`, `sweep_max`, `sweep_steps` — range definition
- All other parameters: fixed values

**View flow:**

1. Generate N values via `numpy.linspace`
2. Enqueue N Celery tasks
3. Wait page; when all complete, render a results grid
4. Plot quantity of interest (e.g., max temperature) vs swept parameter

---

### 5.6 User Accounts and Job History

Use Django's built-in auth system:

```python
from django.contrib.auth.decorators import login_required

@login_required
def conduction(request):
    ...
```

- Add registration / login / logout pages
- Associate every `Mesh`, `SolverJob`, and result with `request.user`
- Users only see their own data
- Django Admin view lets the thesis supervisor review student runs

---

### 5.7 Comprehensive Documentation (JOSS Requirement)

JOSS reviewers check documentation quality explicitly. All items below must exist before Phase 7.

#### README overhaul

Replace the current 4-line README with sections covering:

- **Badges row** — CI status, license, JOSS (added in Phase 7), Python version
- **Screenshot** — actual browser screenshot of the solver in action
- **What it does** — 1–2 sentences on the tool and its physics scope
- **Installation** — conda and Docker paths with copy-pasteable commands
- **Quickstart** — how to run and open the first solver in 5 steps
- **Features table** — solvers, BC types, export formats, visualization options
- **Scientific background** — 1 paragraph on FEM and equations solved; link to `docs/science/governing_equations.md`
- **Validation** — 1 sentence + link to the convergence study page
- **Citation** — BibTeX block pointing to Zenodo DOI (added in Phase 7)
- **License** — MIT

#### API docstrings

Add one-line docstrings to all public functions. Standard format:

```python
def solve_conduction(dimension, mesh, left_bc, right_bc, ...):
    """Solve steady-state heat conduction via FEniCS FEM; return (u_sol, solution_list)."""
```

Functions to cover:

- `SteadyStateThermal/conduction.py` — `solve_conduction`, `generate_solution_plot`
- `ConvectionDiffusion/solver.py` — `solve_convection_diffusion_galerkin`, `solve_convection_diffusion_supg`
- `TransientThermal/solver.py` — `solve_transient`
- `mesh/mesh.py` — `generate_mesh`, `show_mesh`, `generate_interactive_mesh_plot`
- `Validation/convergence.py` — `run_h_refinement_study`, `run_p_refinement_study`

#### `CONTRIBUTING.md` (full upgrade from Phase 0 skeleton)

Expand to include:

- Development environment setup
- How to run tests (`pytest tests/ -v`)
- How to add a new solver (step-by-step, pointing to existing solvers as templates)
- Code style guide (PEP 8, max line length 100)
- How to submit a pull request

#### User guide pages

Add two short pages under `docs/usage/`:

- `getting_started.md` — first solve walkthrough with screenshots
- `solvers.md` — what each solver does, what inputs it expects, what the results mean

---

## Phase 5 Files to Create

| File | Purpose |
| ---- | ------- |
| `Jobs/__init__.py` | App init |
| `Jobs/apps.py` | App config |
| `Jobs/models.py` | SolverJob model |
| `Jobs/tasks.py` | Celery tasks for each solver type |
| `Jobs/views.py` | JobStatusView, JobDashboardView, JobAPIView |
| `Jobs/urls.py` | URL routes |
| `Jobs/templates/jobs/dashboard.html` | Job history table |
| `Jobs/templates/jobs/status.html` | Job polling / waiting page |
| `CONTRIBUTING.md` | Full contributor guide |
| `docs/usage/getting_started.md` | User walkthrough |
| `docs/usage/solvers.md` | Solver reference |

**Files to modify:**

| File | Change |
| ---- | ------ |
| `README.md` | Full overhaul per 5.7 spec |
| `docker-compose.yml` | Add Redis + Celery worker services |
| `Dockerfile` | Ensure Celery is installed |
| `FenicsWeb/settings.py` | Celery config, AUTH, MEDIA_ROOT |
| `FenicsWeb/urls.py` | Jobs URLs, auth URLs |
| `mesh/forms.py` | Uncomment CustomMesh, add FileField |
| All solver `views.py` | Wrap in `@login_required`, enqueue to Celery |
| All solver Python files | Add one-line docstrings to public functions |

---

## Definition of Done

- [ ] Celery + Redis running in Docker, async jobs work end-to-end
- [ ] Job status page polls and redirects to results when done
- [ ] Jobs dashboard shows all user jobs with status
- [ ] VTK export button works, file opens in ParaView
- [ ] XDMF export button works
- [ ] Gmsh `.msh` file upload accepted and converted
- [ ] Parameter sweep page runs N jobs and plots results vs parameter
- [ ] User registration / login / logout working
- [ ] All solver pages require login
- [ ] README has badges, screenshot, installation, quickstart, citation block
- [ ] All public solver functions have one-line docstrings
- [ ] `CONTRIBUTING.md` covers setup, tests, adding a solver, PR process
- [ ] `docs/usage/` has getting-started and solver-reference pages
