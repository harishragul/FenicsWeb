# FenicsWeb — Backend Architecture Plan (UI Redesign Support)

**Companion to:** `docs/ui/app_flow.md`
**Purpose:** Every backend change needed to support the project-based UI — models, refactoring, settings, URLs, solver bridge, file storage.

---

## Current Backend Problems

| Problem | Impact |
| ------- | ------ |
| `generate_mesh(mesh_type, request)` reads `request.POST` directly | Cannot call from project view without a live HTTP request |
| `generate_mesh_plot` / `generate_solution_plot` write to `static/mesh_plot.png` (shared) | Multiple projects overwrite each other's plots |
| All state lives in Django session (`request.session['mesh_str']`) | Lost on session expiry; not queryable; no history |
| `solve_conduction` takes 10+ flat kwargs | Hard to drive from JSON config stored in DB |
| `TEMPLATES['DIRS']` is empty | Shared `base.html` outside any app's templates folder cannot be found |
| Root URL `/` returns 404 | No home page |
| Old `mesh.models.Mesh` is a global pool shared by all users/problems | Not associated with any project |

---

## Architecture After Refactoring

```
HTTP Request
     │
     ▼
projects/views.py  ←── orchestration layer (owns UI, models, routing)
     │
     ├── projects/solver_bridge.py  ←── unpacks ProjectSetup JSON → solver kwargs
     │        │
     │        ▼
     │   SteadyStateThermal/conduction.py   ← pure solver (unchanged API)
     │   ConvectionDiffusion/solver.py      ← Phase 2 (future)
     │   TransientThermal/solver.py         ← Phase 3 (future)
     │
     └── mesh/mesh.py  ←── pure mesh utilities (decoupled from request)
              │
              ▼
         media/projects/<id>/mesh_plot.png
         media/projects/<id>/solution_plot.png
```

**Rule:** Solver libraries (`SteadyStateThermal`, `ConvectionDiffusion`, etc.) never import Django views, models, or `request`. They are pure Python — callable from tests, views, and scripts equally.

---

## Section 1 — New `projects` App

### 1.1 Models (`projects/models.py`)

```python
from django.db import models


class Project(models.Model):
    SOLVER_CHOICES = [
        ('steady_conduction',    'Steady Conduction'),
        ('convection_diffusion', 'Convection-Diffusion'),  # Phase 2
        ('transient',            'Transient Heat'),         # Phase 3
        ('inverse',              'Inverse Problem'),        # Phase 6
    ]
    STATUS = [
        ('new',        'New'),
        ('mesh_done',  'Mesh Complete'),
        ('setup_done', 'Setup Complete'),
        ('solved',     'Solved'),
    ]
    name        = models.CharField(max_length=200)
    solver_type = models.CharField(max_length=50, choices=SOLVER_CHOICES,
                                    default='steady_conduction')
    status      = models.CharField(max_length=20, choices=STATUS, default='new')
    created_at  = models.DateTimeField(auto_now_add=True)
    updated_at  = models.DateTimeField(auto_now=True)

    def __str__(self):
        return f"{self.name} ({self.get_solver_type_display()})"

    def mesh_plot_url(self):
        if hasattr(self, 'mesh_config') and self.mesh_config.mesh_plot_path:
            return f'/media/{self.mesh_config.mesh_plot_path}'
        return None

    def solution_plot_url(self):
        if hasattr(self, 'result') and self.result.solution_plot_path:
            return f'/media/{self.result.solution_plot_path}'
        return None


class ProjectMesh(models.Model):
    project        = models.OneToOneField(Project, on_delete=models.CASCADE,
                                           related_name='mesh_config')
    mesh_type      = models.CharField(max_length=50)
    mesh_str       = models.TextField()        # serialised for show_mesh()
    mesh_params    = models.JSONField()        # human-readable param summary
    mesh_plot_path = models.CharField(max_length=300, blank=True)
                                               # relative to MEDIA_ROOT


class ProjectSetup(models.Model):
    project       = models.OneToOneField(Project, on_delete=models.CASCADE,
                                          related_name='setup_config')
    dimension     = models.CharField(max_length=5)   # '1D', '2D', '3D'
    bc_config     = models.JSONField()               # per-face BC type + values
    solver_params = models.JSONField()               # k, f, element_order, …


class ProjectResult(models.Model):
    project             = models.OneToOneField(Project, on_delete=models.CASCADE,
                                               related_name='result')
    mesh_plot_path      = models.CharField(max_length=300)
    solution_plot_path  = models.CharField(max_length=300)
    solution_data       = models.JSONField()   # vertex coords + temperature values
    solved_at           = models.DateTimeField(auto_now_add=True)
    solve_time_s        = models.FloatField(null=True)
```

**`bc_config` JSON schema** (stored in `ProjectSetup`):

```json
{
  "left":   { "type": "dirichlet", "value": 0.0 },
  "right":  { "type": "robin",     "h": 10.0, "u_inf": 20.0 },
  "bottom": { "type": "dirichlet", "value": 0.0 },
  "top":    { "type": "neumann",   "value": 5.0 }
}
```

**`solver_params` JSON schema:**

```json
{ "k": 1.0, "f": 0.0, "element_order": 1 }
```

---

### 1.2 Solver Bridge (`projects/solver_bridge.py`)

This is the key new backend piece. It unpacks the project JSON config and calls the appropriate solver library function.

```python
import time
from mesh.mesh import show_mesh
from SteadyStateThermal.conduction import solve_conduction


def run_project_solver(project):
    """Run the solver for a project; return (u_sol, solution_list, solve_time_s)."""
    setup = project.setup_config
    mesh_cfg = project.mesh_config

    mesh = show_mesh(mesh_cfg.mesh_str)
    bc   = setup.bc_config
    sp   = setup.solver_params

    if project.solver_type == 'steady_conduction':
        return _run_steady_conduction(mesh, setup.dimension, bc, sp)

    # Future solvers plugged in here (Phase 2, 3, 6):
    # elif project.solver_type == 'convection_diffusion':
    #     from ConvectionDiffusion.solver import solve_convection_diffusion_supg
    #     return _run_convection_diffusion(mesh, setup.dimension, bc, sp)

    raise ValueError(f"Unknown solver type: {project.solver_type}")


def _run_steady_conduction(mesh, dimension, bc, sp):
    bc_types      = {face: cfg['type']         for face, cfg in bc.items()}
    neumann_fluxes = {face: cfg.get('value', 0) for face, cfg in bc.items()
                      if cfg['type'] == 'neumann'}
    robin_hs      = {face: cfg.get('h', 0)     for face, cfg in bc.items()
                      if cfg['type'] == 'robin'}
    robin_u_infs  = {face: cfg.get('u_inf', 0) for face, cfg in bc.items()
                      if cfg['type'] == 'robin'}

    def _v(face):
        cfg = bc.get(face, {})
        return cfg.get('value', 0.0) if cfg.get('type') == 'dirichlet' else 0.0

    t0 = time.time()
    u_sol, solution_list = solve_conduction(
        dimension, mesh,
        left_bc=_v('left'), right_bc=_v('right'),
        top_bc=_v('top'),   bottom_bc=_v('bottom'),
        front_bc=_v('front'), back_bc=_v('back'),
        f=sp.get('f', 0.0),
        k=sp.get('k', 1.0),
        bc_types=bc_types,
        neumann_fluxes=neumann_fluxes,
        robin_hs=robin_hs,
        robin_u_infs=robin_u_infs,
    )
    return u_sol, solution_list, time.time() - t0
```

---

### 1.3 URL Routes (`projects/urls.py`)

```python
from django.urls import path
from . import views

urlpatterns = [
    path('',                    views.home,       name='home'),
    path('new/',                views.create,     name='project_create'),
    path('<int:pk>/',           views.overview,   name='project_overview'),
    path('<int:pk>/mesh/',      views.mesh_step,  name='project_mesh'),
    path('<int:pk>/setup/',     views.setup_step, name='project_setup'),
    path('<int:pk>/solve/',     views.solve_step, name='project_solve'),
    path('<int:pk>/delete/',    views.delete,     name='project_delete'),
]
```

---

### 1.4 Forms (`projects/forms.py`)

```python
from django import forms
from .models import Project

MESH_TYPE_CHOICES = [
    ('interval',      'IntervalMesh (1D)'),
    ('rectangle',     'RectangleMesh (2D)'),
    ('box',           'BoxMesh (3D)'),
    ('unit_interval', 'UnitIntervalMesh (1D)'),
    ('unit_square',   'UnitSquareMesh (2D)'),
    ('unit_cube',     'UnitCubeMesh (3D)'),
]

class ProjectCreateForm(forms.ModelForm):
    class Meta:
        model  = Project
        fields = ['name', 'solver_type']
        widgets = {
            'name':        forms.TextInput(attrs={'class': 'form-control', 'placeholder': 'e.g. Rod Analysis'}),
            'solver_type': forms.Select(attrs={'class': 'form-select'}),
        }

class MeshStepForm(forms.Form):
    mesh_type = forms.ChoiceField(choices=MESH_TYPE_CHOICES,
                                   widget=forms.Select(attrs={'class': 'form-select', 'id': 'mesh_type_select'}))
    # Dynamic geometry fields — all optional; view validates based on mesh_type
    n       = forms.IntegerField(required=False, label='Number of cells (n)')
    x0      = forms.FloatField(required=False,   label='x₀')
    x1      = forms.FloatField(required=False,   label='x₁')
    y0      = forms.FloatField(required=False,   label='y₀')
    y1      = forms.FloatField(required=False,   label='y₁')
    z0      = forms.FloatField(required=False,   label='z₀')
    z1      = forms.FloatField(required=False,   label='z₁')
    nx      = forms.IntegerField(required=False, label='nx')
    ny      = forms.IntegerField(required=False, label='ny')
    nz      = forms.IntegerField(required=False, label='nz')

class SetupStepForm(forms.Form):
    BC_TYPES = [
        ('dirichlet', 'Dirichlet (fixed T)'),
        ('neumann',   'Neumann (heat flux)'),
        ('robin',     'Robin (convective cooling)'),
    ]
    DIMS = [('1D', '1D'), ('2D', '2D'), ('3D', '3D')]

    dimension = forms.ChoiceField(choices=DIMS, widget=forms.Select(attrs={'class': 'form-select'}))

    # One set of fields per face; view renders only the active set based on dimension
    left_type   = forms.ChoiceField(choices=BC_TYPES)
    left_value  = forms.FloatField(required=False, initial=0.0)
    left_h      = forms.FloatField(required=False)
    left_u_inf  = forms.FloatField(required=False)

    right_type  = forms.ChoiceField(choices=BC_TYPES)
    right_value = forms.FloatField(required=False, initial=0.0)
    right_h     = forms.FloatField(required=False)
    right_u_inf = forms.FloatField(required=False)

    bottom_type  = forms.ChoiceField(choices=BC_TYPES)
    bottom_value = forms.FloatField(required=False, initial=0.0)
    bottom_h     = forms.FloatField(required=False)
    bottom_u_inf = forms.FloatField(required=False)

    top_type  = forms.ChoiceField(choices=BC_TYPES)
    top_value = forms.FloatField(required=False, initial=0.0)
    top_h     = forms.FloatField(required=False)
    top_u_inf = forms.FloatField(required=False)

    front_type  = forms.ChoiceField(choices=BC_TYPES)
    front_value = forms.FloatField(required=False, initial=0.0)
    front_h     = forms.FloatField(required=False)
    front_u_inf = forms.FloatField(required=False)

    back_type  = forms.ChoiceField(choices=BC_TYPES)
    back_value = forms.FloatField(required=False, initial=0.0)
    back_h     = forms.FloatField(required=False)
    back_u_inf = forms.FloatField(required=False)

    k = forms.FloatField(label='Thermal conductivity k (W/m·K)', initial=1.0)
    f = forms.FloatField(label='Source term f',                   initial=0.0)
```

---

### 1.5 Views (`projects/views.py`) — logic outline

```python
# home — list all projects
def home(request):
    projects = Project.objects.all().order_by('-updated_at')
    return render(request, 'projects/home.html', {'projects': projects})

# create — new project form
def create(request):
    if request.method == 'POST':
        form = ProjectCreateForm(request.POST)
        if form.is_valid():
            project = form.save()
            return redirect('project_mesh', pk=project.pk)
    else:
        form = ProjectCreateForm()
    return render(request, 'projects/create.html', {'form': form})

# overview — pipeline view
def overview(request, pk):
    project = get_object_or_404(Project, pk=pk)
    return render(request, 'projects/overview.html', {'project': project})

# mesh_step — define mesh, save ProjectMesh, generate preview
def mesh_step(request, pk):
    project = get_object_or_404(Project, pk=pk)
    if request.method == 'POST':
        form = MeshStepForm(request.POST)
        if form.is_valid():
            mesh, mesh_str, mesh_params = _build_mesh_from_form(form)
            output_dir = _project_media_dir(project)
            plot_rel   = generate_mesh_plot(mesh, output_dir=output_dir)
            ProjectMesh.objects.update_or_create(
                project=project,
                defaults={
                    'mesh_type':      form.cleaned_data['mesh_type'],
                    'mesh_str':       mesh_str,
                    'mesh_params':    mesh_params,
                    'mesh_plot_path': plot_rel,
                },
            )
            project.status = 'mesh_done'
            project.save()
            return redirect('project_setup', pk=project.pk)
    else:
        initial = {}
        if hasattr(project, 'mesh_config'):
            initial = project.mesh_config.mesh_params
        form = MeshStepForm(initial=initial)
    return render(request, 'projects/mesh.html', {'project': project, 'form': form})

# setup_step — define BCs, save ProjectSetup
def setup_step(request, pk):
    project = get_object_or_404(Project, pk=pk)
    if request.method == 'POST':
        form = SetupStepForm(request.POST)
        if form.is_valid():
            bc_config, solver_params, dimension = _extract_setup(form)
            ProjectSetup.objects.update_or_create(
                project=project,
                defaults={
                    'dimension':     dimension,
                    'bc_config':     bc_config,
                    'solver_params': solver_params,
                },
            )
            project.status = 'setup_done'
            project.save()
            return redirect('project_overview', pk=project.pk)
    else:
        initial = {}
        if hasattr(project, 'setup_config'):
            initial = _setup_to_form_initial(project.setup_config)
        form = SetupStepForm(initial=initial)
    return render(request, 'projects/setup.html', {'project': project, 'form': form})

# solve_step — run solver, save ProjectResult
def solve_step(request, pk):
    project = get_object_or_404(Project, pk=pk)
    if request.method == 'POST':
        u_sol, solution_list, elapsed = run_project_solver(project)
        output_dir    = _project_media_dir(project)
        mesh          = show_mesh(project.mesh_config.mesh_str)
        sol_plot_rel  = generate_solution_plot(u_sol, mesh, output_dir=output_dir)
        mesh_plot_rel = generate_mesh_plot(mesh, output_dir=output_dir)
        ProjectResult.objects.update_or_create(
            project=project,
            defaults={
                'mesh_plot_path':     mesh_plot_rel,
                'solution_plot_path': sol_plot_rel,
                'solution_data':      solution_list,
                'solve_time_s':       elapsed,
            },
        )
        project.status = 'solved'
        project.save()
        return redirect('project_solve', pk=project.pk)
    return render(request, 'projects/solve.html', {'project': project})

# delete
def delete(request, pk):
    project = get_object_or_404(Project, pk=pk)
    if request.method == 'POST':
        project.delete()
        return redirect('home')
    return render(request, 'projects/confirm_delete.html', {'project': project})


# ── Helpers ──────────────────────────────────────────────────────────────────

def _project_media_dir(project):
    """Return (and create) per-project media subdirectory path."""
    from django.conf import settings
    path = os.path.join(settings.MEDIA_ROOT, 'projects', str(project.pk))
    os.makedirs(path, exist_ok=True)
    return path

def _build_mesh_from_form(form):
    """Extract validated form data → (FEniCS mesh, mesh_str, params_dict)."""
    d = form.cleaned_data
    mt = d['mesh_type']
    if mt == 'interval':
        mesh = IntervalMesh(d['n'], d['x0'], d['x1'])
        mesh_str = f"interval,{d['n']},{d['x0']},{d['x1']}"
        params = {'mesh_type': mt, 'n': d['n'], 'x0': d['x0'], 'x1': d['x1']}
    elif mt == 'unit_interval':
        mesh = UnitIntervalMesh(d['n'])
        mesh_str = f"unit_interval,{d['n']}"
        params = {'mesh_type': mt, 'n': d['n']}
    elif mt == 'rectangle':
        mesh = RectangleMesh(Point(d['x0'], d['y0']), Point(d['x1'], d['y1']), d['nx'], d['ny'])
        mesh_str = f"rectangle,{d['x0']},{d['y0']},{d['x1']},{d['y1']},{d['nx']},{d['ny']}"
        params = {'mesh_type': mt, 'x0': d['x0'], 'y0': d['y0'],
                  'x1': d['x1'], 'y1': d['y1'], 'nx': d['nx'], 'ny': d['ny']}
    # … (box, unit_square, unit_cube similarly)
    return mesh, mesh_str, params

def _extract_setup(form):
    """Convert SetupStepForm → (bc_config dict, solver_params dict, dimension)."""
    d = form.cleaned_data
    faces = ['left', 'right', 'bottom', 'top', 'front', 'back']
    bc_config = {}
    for face in faces:
        ftype = d.get(f'{face}_type', 'dirichlet')
        entry = {'type': ftype}
        if ftype == 'dirichlet':
            entry['value'] = d.get(f'{face}_value') or 0.0
        elif ftype == 'neumann':
            entry['value'] = d.get(f'{face}_value') or 0.0
        elif ftype == 'robin':
            entry['h']     = d.get(f'{face}_h')     or 0.0
            entry['u_inf'] = d.get(f'{face}_u_inf') or 0.0
        bc_config[face] = entry
    solver_params = {'k': d['k'], 'f': d['f']}
    return bc_config, solver_params, d['dimension']
```

---

## Section 2 — Solver Library Refactoring

### 2.1 `mesh/mesh.py` — decouple from `request`

**Problem:** `generate_mesh(mesh_type, request)` reads `request.POST` directly, making it unusable outside a Django view context.

**Fix:** Add `build_mesh_from_params(mesh_type, params)` — a pure function that takes a dict. The existing `generate_mesh(mesh_type, request)` calls it internally for backward compatibility.

```python
def build_mesh_from_params(mesh_type, params):
    """Create FEniCS mesh from a plain dict — no Django request needed."""
    if mesh_type == 'interval':
        return (IntervalMesh(params['n'], params['x0'], params['x1']),
                f"interval,{params['n']},{params['x0']},{params['x1']}")
    elif mesh_type == 'unit_interval':
        return UnitIntervalMesh(params['n']), f"unit_interval,{params['n']}"
    elif mesh_type == 'rectangle':
        return (RectangleMesh(Point(params['x0'], params['y0']),
                              Point(params['x1'], params['y1']),
                              params['nx'], params['ny']),
                f"rectangle,{params['x0']},{params['y0']},"
                f"{params['x1']},{params['y1']},{params['nx']},{params['ny']}")
    # … box, unit_square, unit_cube
    raise ValueError(f"Unknown mesh type: {mesh_type}")

def generate_mesh(mesh_type, request):
    """Existing entry point — extract params from request.POST, call build_mesh_from_params."""
    # (current form-parsing logic stays here unchanged)
    ...
```

**Backward-compatibility:** existing `generate_mesh(mesh_type, request)` keeps working unchanged. Old tests pass.

---

### 2.2 `mesh/mesh.py` and `conduction.py` — per-project output paths

Both plotting functions currently hardcode their output path to `static/`. Add an `output_dir` parameter defaulting to `None` (which falls back to the current `STATIC_DIR`).

**`generate_mesh_plot` signature change:**

```python
def generate_mesh_plot(mesh, output_dir=None):
    if output_dir is None:
        output_dir = STATIC_DIR
    os.makedirs(output_dir, exist_ok=True)
    plot_filename = os.path.join(output_dir, 'mesh_plot.png')
    # … rest unchanged
    return 'mesh_plot.png'   # relative filename only
```

**`generate_solution_plot` signature change:**

```python
def generate_solution_plot(u_sol, mesh, output_dir=None):
    if output_dir is None:
        output_dir = STATIC_DIR
    os.makedirs(output_dir, exist_ok=True)
    plot_filename = os.path.join(output_dir, 'solution_plot.png')
    # … rest unchanged
    return 'solution_plot.png'
```

**Backward-compatibility:** all existing calls `generate_mesh_plot(mesh)` and `generate_solution_plot(u_sol, mesh)` continue to write to `static/` as before. All existing tests pass with no changes.

---

## Section 3 — Django Settings Changes (`FenicsWeb/settings.py`)

### 3.1 `INSTALLED_APPS`

```python
INSTALLED_APPS = [
    ...
    'main',
    'mesh',
    'SteadyStateThermal',
    'projects',          # ADD
]
```

### 3.2 `TEMPLATES['DIRS']`

```python
from pathlib import Path
BASE_DIR = Path(__file__).resolve().parent.parent

TEMPLATES = [
    {
        ...
        'DIRS': [BASE_DIR / 'templates'],   # ADD — enables shared base.html
        'APP_DIRS': True,
        ...
    },
]
```

### 3.3 Media files (for per-project plots)

```python
MEDIA_URL  = '/media/'
MEDIA_ROOT = BASE_DIR / 'media'
```

Generated plots go to `media/projects/<id>/mesh_plot.png` instead of `static/`. This is the correct Django pattern — `static/` is for committed assets, `media/` is for generated/user content.

---

## Section 4 — URL Structure Changes (`FenicsWeb/urls.py`)

```python
from django.contrib import admin
from django.urls import path, include
from django.conf import settings
from django.conf.urls.static import static

urlpatterns = [
    path('admin/',                admin.site.urls),
    path('',                      include('projects.urls')),            # NEW — root = home
    path('SteadyStateThermal/',   include('SteadyStateThermal.urls')), # keep (backward compat)
    path('mesh/',                 include('mesh.urls')),                # keep (backward compat)
    path('validation/',           include('SteadyStateThermal.urls')), # keep
] + static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)      # NEW — serve media
```

---

## Section 5 — Shared Template Structure

```
templates/                          ← BASE_DIR/templates/ (new top-level dir)
├── base.html                       ← sidebar + navbar layout
├── partials/
│   ├── sidebar.html                ← solver list with active state
│   ├── pipeline_bar.html           ← Mesh → Setup → Solve indicator
│   └── project_card.html           ← reusable project summary card
└── projects/
    ├── home.html
    ├── create.html
    ├── overview.html
    ├── mesh.html
    ├── setup.html
    ├── solve.html
    └── confirm_delete.html
```

`base.html` skeleton:

```html
<!DOCTYPE html>
<html>
<head>
    <title>FenicsWeb {% block page_title %}{% endblock %}</title>
    <!-- Bootstrap 5, Plotly CDN links -->
</head>
<body>
<div class="d-flex" style="min-height:100vh">

    <!-- Sidebar -->
    <div class="bg-dark text-white" style="width:240px; min-height:100vh">
        {% include 'partials/sidebar.html' %}
    </div>

    <!-- Main area -->
    <div class="flex-grow-1 d-flex flex-column">
        <!-- Navbar -->
        <nav class="navbar navbar-light bg-light px-3">
            {% block breadcrumb %}{% endblock %}
        </nav>
        <!-- Content -->
        <div class="p-4 flex-grow-1">
            {% block content %}{% endblock %}
        </div>
    </div>
</div>
</body>
</html>
```

---

## Section 6 — Test Updates

### 6.1 New tests for `projects` app (`tests/test_projects.py`)

```python
import pytest
from django.test import TestCase, Client
from projects.models import Project, ProjectMesh, ProjectSetup


class TestProjectModels(TestCase):
    def test_create_project(self):
        p = Project.objects.create(name='Test Rod', solver_type='steady_conduction')
        assert p.status == 'new'
        assert str(p) == 'Test Rod (Steady Conduction)'

    def test_project_status_progression(self):
        p = Project.objects.create(name='Test', solver_type='steady_conduction')
        assert p.status == 'new'
        p.status = 'mesh_done'
        p.save()
        assert Project.objects.get(pk=p.pk).status == 'mesh_done'


class TestSolverBridge(TestCase):
    def test_bridge_steady_conduction_1d(self):
        from projects.solver_bridge import run_project_solver
        p = Project.objects.create(name='Bridge Test', solver_type='steady_conduction')
        ProjectMesh.objects.create(
            project=p,
            mesh_type='unit_interval',
            mesh_str='unit_interval,10',
            mesh_params={'n': 10},
        )
        ProjectSetup.objects.create(
            project=p,
            dimension='1D',
            bc_config={
                'left':  {'type': 'dirichlet', 'value': 0.0},
                'right': {'type': 'dirichlet', 'value': 100.0},
            },
            solver_params={'k': 1.0, 'f': 0.0},
        )
        u_sol, solution_list, elapsed = run_project_solver(p)
        assert u_sol is not None
        assert len(solution_list) > 0
        assert elapsed > 0
        # Check midpoint temperature ≈ 50 (linear interpolation)
        mid_temp = next(t for (coords, t) in solution_list if abs(coords[0] - 0.5) < 0.06)
        assert abs(mid_temp - 50.0) < 1.0


class TestProjectViews(TestCase):
    def setUp(self):
        self.client = Client()

    def test_home_page_loads(self):
        response = self.client.get('/')
        assert response.status_code == 200

    def test_create_project_redirects_to_mesh(self):
        response = self.client.post('/new/', {
            'name': 'Test', 'solver_type': 'steady_conduction'
        })
        assert response.status_code == 302
        assert '/mesh/' in response['Location']
```

### 6.2 Updated `test_conduction.py` — test `output_dir` parameter

```python
def test_generate_solution_plot_writes_to_custom_dir(tmp_path):
    from SteadyStateThermal.conduction import solve_conduction, generate_solution_plot
    from fenics import UnitIntervalMesh
    mesh = UnitIntervalMesh(10)
    u_sol, _ = solve_conduction('1D', mesh, left_bc=0, right_bc=1,
                                top_bc=0, bottom_bc=0, front_bc=0, back_bc=0,
                                f=0, k=1)
    rel = generate_solution_plot(u_sol, mesh, output_dir=str(tmp_path))
    assert (tmp_path / rel).exists()

def test_generate_mesh_plot_writes_to_custom_dir(tmp_path):
    from mesh.mesh import generate_mesh_plot
    from fenics import UnitSquareMesh
    mesh = UnitSquareMesh(4, 4)
    rel = generate_mesh_plot(mesh, output_dir=str(tmp_path))
    assert (tmp_path / rel).exists()
```

---

## Section 7 — Migration Plan (No Breaking Changes)

All changes preserve backward compatibility. Old URLs and tests keep working throughout.

| Step | Action | Breaking? |
| ---- | ------ | --------- |
| 1 | Create `projects` app, models, migrations | No |
| 2 | Add `output_dir` param to plot functions (default=None) | No — old calls unchanged |
| 3 | Add `build_mesh_from_params` to `mesh/mesh.py` | No — additive |
| 4 | Add `MEDIA_ROOT`, `MEDIA_URL`, `TEMPLATES DIRS` to settings | No |
| 5 | Add `projects.urls` to root URL conf | No — old URLs still work |
| 6 | Build templates (`base.html`, sidebar, partials) | No |
| 7 | Build home/create/overview/mesh/setup/solve views | No |
| 8 | Wire navigation links (navbar, sidebar) | No |
| 9 | Deprecate old `Mesh` model (keep table, stop writing to it) | No |
| 10 | Update `tests/` with project + bridge tests | No |

---

## Complete File Change List

### New files

| File | Purpose |
| ---- | ------- |
| `projects/__init__.py` | App init |
| `projects/apps.py` | AppConfig |
| `projects/models.py` | Project, ProjectMesh, ProjectSetup, ProjectResult |
| `projects/views.py` | home, create, overview, mesh_step, setup_step, solve_step, delete |
| `projects/forms.py` | ProjectCreateForm, MeshStepForm, SetupStepForm |
| `projects/urls.py` | URL routing |
| `projects/solver_bridge.py` | Unpack project config → call solver |
| `projects/migrations/0001_initial.py` | Auto-generated migration |
| `templates/base.html` | Shared sidebar + navbar layout |
| `templates/partials/sidebar.html` | Solver list |
| `templates/partials/pipeline_bar.html` | Mesh → Setup → Solve indicator |
| `templates/partials/project_card.html` | Project summary card |
| `templates/projects/home.html` | Main page |
| `templates/projects/create.html` | New project form |
| `templates/projects/overview.html` | Pipeline view |
| `templates/projects/mesh.html` | Mesh step |
| `templates/projects/setup.html` | Setup step |
| `templates/projects/solve.html` | Results page |
| `templates/projects/confirm_delete.html` | Delete confirmation |
| `tests/test_projects.py` | Project model, bridge, and view tests |

### Modified files

| File | Change |
| ---- | ------ |
| `FenicsWeb/settings.py` | Add `projects` to INSTALLED_APPS, TEMPLATES DIRS, MEDIA_ROOT/URL |
| `FenicsWeb/urls.py` | Add projects root URL, media URL serving |
| `mesh/mesh.py` | Add `build_mesh_from_params()`; add `output_dir` to `generate_mesh_plot()` |
| `SteadyStateThermal/conduction.py` | Add `output_dir` to `generate_solution_plot()` |
| `tests/test_conduction.py` | Add `output_dir` tests using `tmp_path` fixture |
