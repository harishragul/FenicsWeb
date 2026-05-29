# FenicsWeb — UI Architecture & App Flow Plan

**Why this exists before Phase 2:**
The current app uses Django sessions to carry state between pages, has no concept of a "project,"
and each solver lives in its own disconnected URL space. Adding more solvers on top of this
structure would create an unmaintainable tangle of session keys and URLs.
This document defines the new project-based architecture that all future phases build on.

---

## Target Experience (ANSYS Workbench Model)

```
┌─────────────────┬──────────────────────────────────────────────────┐
│  SIDEBAR        │  Main Area                                        │
│                 │                                                    │
│  ● Solvers      │  My Projects                        [+ New]       │
│                 │                                                    │
│  ○ Steady       │  ┌──────────────────────────────────────────────┐ │
│    Conduction   │  │  Rod Analysis          Steady Conduction     │ │
│                 │  │  [Mesh ✓]──[Setup ✓]──[Solve ○]   [Open ▶] │ │
│  ○ Convection-  │  └──────────────────────────────────────────────┘ │
│    Diffusion    │                                                    │
│    (Phase 2)    │  ┌──────────────────────────────────────────────┐ │
│                 │  │  Heat Sink 2D       Convection-Diffusion     │ │
│  ○ Transient    │  │  [Mesh ✓]──[Setup ○]──[Solve ○]   [Open ▶] │ │
│    (Phase 3)    │  └──────────────────────────────────────────────┘ │
│                 │                                                    │
│  ── Tools ──    │                                                    │
│  ○ Validation   │                                                    │
│  ○ Inverse      │                                                    │
│    (Phase 6)    │                                                    │
└─────────────────┴──────────────────────────────────────────────────┘
```

---

## App Flow Diagram

```
[Main Page /]
     │
     ├── [+ New Project] → Choose solver type → Name project → [Project Overview]
     │
     └── [Open ▶ existing project] ──────────────────────────────────────┐
                                                                          ▼
                                                              [Project Overview /project/<id>/]
                                                              Shows pipeline bar:
                                                              Mesh → Setup → Solve
                                                                    │
                                          ┌───────────────────────┼───────────────────────┐
                                          ▼                       ▼                       ▼
                              [Mesh Step]              [Setup Step]            [Solve Step]
                              /project/<id>/mesh/      /project/<id>/setup/    /project/<id>/solve/
                              Define geometry          Define BCs + params     Run + view results
                              Save → back to           Save → back to          Download plots/CSV
                              overview                 overview
```

---

## Screen Layouts

### Screen 1 — Main Page (`/`)

```
┌────────────────────────────────────────────────────────────────────┐
│  🔥 FenicsWeb                                              [About] │
├─────────────────┬──────────────────────────────────────────────────┤
│                 │                                                   │
│  SOLVERS        │  Projects                            [+ New]     │
│  ─────────────  │                                                   │
│  ● Steady       │  ╔══════════════════════════════════════════════╗ │
│    Conduction   │  ║  Rod Analysis                               ║ │
│                 │  ║  Solver: Steady Conduction                  ║ │
│  ○ Convection-  │  ║                                             ║ │
│    Diffusion    │  ║  [Mesh ✓]──[Setup ✓]──[Solve ✓]            ║ │
│    (soon)       │  ║  Last solved: 2026-05-29                    ║ │
│                 │  ║                                    [Open ▶] ║ │
│  ○ Transient    │  ╚══════════════════════════════════════════════╝ │
│    (soon)       │                                                   │
│                 │  ╔══════════════════════════════════════════════╗ │
│  ─────────────  │  ║  Heat Sink Study                            ║ │
│  TOOLS          │  ║  Solver: Steady Conduction                  ║ │
│                 │  ║                                             ║ │
│  ○ Validation   │  ║  [Mesh ✓]──[Setup ○]──[Solve ○]            ║ │
│  ○ Inverse      │  ║  Setup not complete                         ║ │
│    Problem      │  ║                                    [Open ▶] ║ │
│    (Phase 6)    │  ╚══════════════════════════════════════════════╝ │
│                 │                                                   │
└─────────────────┴───────────────────────────────────────────────────┘
```

### Screen 2 — Project Overview (`/project/<id>/`)

```
┌────────────────────────────────────────────────────────────────────┐
│  🔥 FenicsWeb                                                      │
├─────────────────┬──────────────────────────────────────────────────┤
│  [← Projects]   │  Rod Analysis                  Steady Conduction │
│                 │  ────────────────────────────────────────────── │
│  SOLVERS        │                                                   │
│  ─────────────  │  Pipeline                                        │
│  ...            │                                                   │
│                 │  ┌──────────┐    ┌──────────┐    ┌──────────┐   │
│                 │  │  MESH    │───▶│  SETUP   │───▶│  SOLVE   │   │
│                 │  │    ✓     │    │    ✓     │    │    ○     │   │
│                 │  │ [Edit]   │    │ [Edit]   │    │ [Run ▶]  │   │
│                 │  └──────────┘    └──────────┘    └──────────┘   │
│                 │                                                   │
│                 │  ── Mesh Summary ──────────────────────────────  │
│                 │  Type: IntervalMesh (1D)                         │
│                 │  n=20, domain [0, 1]                             │
│                 │                                                   │
│                 │  ── Setup Summary ─────────────────────────────  │
│                 │  Left: Dirichlet  T = 0 K                        │
│                 │  Right: Dirichlet  T = 100 K                     │
│                 │  k = 1.0 W/m·K    f = 0                         │
│                 │                                                   │
│                 │  ── Results ───────────────────────────────────  │
│                 │  Not yet solved. Click [Run ▶] to solve.         │
└─────────────────┴───────────────────────────────────────────────────┘
```

### Screen 3 — Mesh Step (`/project/<id>/mesh/`)

```
┌────────────────────────────────────────────────────────────────────┐
│  🔥 FenicsWeb                                                      │
├─────────────────┬──────────────────────────────────────────────────┤
│  [← Overview]   │  Rod Analysis  ▶  Mesh                          │
│                 │  ────────────────────────────────────────────── │
│  SOLVERS        │                                                   │
│  ...            │  Mesh Type ▼                                     │
│                 │  ● IntervalMesh (1D)                             │
│                 │  ○ RectangleMesh (2D)                            │
│                 │  ○ BoxMesh (3D)                                  │
│                 │  ○ Upload Gmsh (.msh)                            │
│                 │                                                   │
│                 │  ┌─────────────────────────────────────────────┐ │
│                 │  │  Number of intervals (n)    [  20  ]       │ │
│                 │  │  Start point x0             [  0   ]       │ │
│                 │  │  End point x1               [  1   ]       │ │
│                 │  └─────────────────────────────────────────────┘ │
│                 │                                                   │
│                 │  [Preview Mesh]              [Save & Continue ▶] │
│                 │                                                   │
│                 │  ┌─────────────────────────────────────────────┐ │
│                 │  │  [mesh preview image renders here]          │ │
│                 │  └─────────────────────────────────────────────┘ │
└─────────────────┴───────────────────────────────────────────────────┘
```

### Screen 4 — Setup Step (`/project/<id>/setup/`)

```
┌────────────────────────────────────────────────────────────────────┐
│  🔥 FenicsWeb                                                      │
├─────────────────┬──────────────────────────────────────────────────┤
│  [← Overview]   │  Rod Analysis  ▶  Setup                         │
│                 │  ────────────────────────────────────────────── │
│  SOLVERS        │                                                   │
│  ...            │  Boundary Conditions                             │
│                 │                                                   │
│                 │  Left face    [Dirichlet ▼]  Value [ 0   ] K    │
│                 │  Right face   [Robin     ▼]  h [10] T∞ [20] K  │
│                 │                                                   │
│                 │  Solver Parameters                               │
│                 │  Thermal conductivity k  [ 1.0  ] W/m·K         │
│                 │  Source term f           [ 0.0  ]               │
│                 │                                                   │
│                 │  [← Back to Mesh]           [Save & Continue ▶] │
└─────────────────┴───────────────────────────────────────────────────┘
```

### Screen 5 — Solve Step (`/project/<id>/solve/`)

```
┌────────────────────────────────────────────────────────────────────┐
│  🔥 FenicsWeb                                                      │
├─────────────────┬──────────────────────────────────────────────────┤
│  [← Overview]   │  Rod Analysis  ▶  Results                       │
│                 │  ────────────────────────────────────────────── │
│  SOLVERS        │                                                   │
│  ...            │  ┌─────────────────┐  ┌────────────────────────┐ │
│                 │  │  Mesh           │  │  Temperature Field     │ │
│                 │  │  [mesh_plot]    │  │  [solution_plot]       │ │
│                 │  └─────────────────┘  └────────────────────────┘ │
│                 │                                                   │
│                 │  [Download CSV]  [Download VTK]  [Re-Run ▶]      │
│                 │                                                   │
│                 │  Solution Data  [Search...]                      │
│                 │  ┌─────────────────────────────────────────────┐ │
│                 │  │  x         Temperature (K)                  │ │
│                 │  │  0.00      0.00                             │ │
│                 │  │  0.05      5.00                             │ │
│                 │  │  ...       ...                              │ │
│                 │  └─────────────────────────────────────────────┘ │
└─────────────────┴───────────────────────────────────────────────────┘
```

---

## Data Models

### New `projects` app

```python
class Project(models.Model):
    SOLVER_CHOICES = [
        ('steady_conduction',     'Steady Conduction'),
        ('convection_diffusion',  'Convection-Diffusion'),   # Phase 2
        ('transient',             'Transient Heat'),          # Phase 3
        ('inverse',               'Inverse Problem'),         # Phase 6
    ]
    STATUS = [
        ('new',        'New'),
        ('mesh_done',  'Mesh Complete'),
        ('setup_done', 'Setup Complete'),
        ('solved',     'Solved'),
    ]
    name        = models.CharField(max_length=200)
    solver_type = models.CharField(max_length=50, choices=SOLVER_CHOICES)
    status      = models.CharField(max_length=20, choices=STATUS, default='new')
    created_at  = models.DateTimeField(auto_now_add=True)
    updated_at  = models.DateTimeField(auto_now=True)


class ProjectMesh(models.Model):
    project     = models.OneToOneField(Project, on_delete=models.CASCADE,
                                       related_name='mesh_config')
    mesh_type   = models.CharField(max_length=50)
    mesh_str    = models.TextField()         # serialised for show_mesh()
    mesh_params = models.JSONField()         # human-readable param summary


class ProjectSetup(models.Model):
    project     = models.OneToOneField(Project, on_delete=models.CASCADE,
                                       related_name='setup_config')
    dimension   = models.CharField(max_length=5)   # '1D', '2D', '3D'
    bc_config   = models.JSONField()               # bc_types, values, h, u_inf per face
    solver_params = models.JSONField()             # k, f, dt, T, scheme, …


class ProjectResult(models.Model):
    project        = models.OneToOneField(Project, on_delete=models.CASCADE,
                                          related_name='result')
    mesh_plot      = models.CharField(max_length=300)
    solution_plot  = models.CharField(max_length=300)
    solution_data  = models.JSONField()    # vertex coords + temperature values
    solved_at      = models.DateTimeField(auto_now_add=True)
    solve_time_s   = models.FloatField(null=True)
```

---

## URL Structure

```
/                              → projects.views.home       (main page)
/project/new/                  → projects.views.create
/project/<id>/                 → projects.views.overview
/project/<id>/mesh/            → projects.views.mesh_step
/project/<id>/setup/           → projects.views.setup_step
/project/<id>/solve/           → projects.views.solve_step
/project/<id>/delete/          → projects.views.delete
/validation/                   → SteadyStateThermal.views.validation  (keep)
```

The old URLs `/SteadyStateThermal/` and `/mesh/` remain accessible during
transition but eventually redirect to the project-based flow.

---

## Template Architecture

```
templates/
├── base.html            ← shared layout: navbar + sidebar + content area
├── projects/
│   ├── home.html        ← project cards grid
│   ├── create.html      ← new project form (name + solver type)
│   ├── overview.html    ← pipeline bar + step summaries
│   ├── mesh.html        ← mesh form + preview
│   ├── setup.html       ← BC form (type-aware per solver)
│   └── solve.html       ← plots + data table + download buttons
└── partials/
    ├── sidebar.html     ← solver list (rendered by base.html)
    ├── pipeline_bar.html ← Mesh → Setup → Solve progress indicator
    └── project_card.html ← reusable project summary card
```

`base.html` includes the sidebar and a `{% block content %}` slot. Every page
extends `base.html`. The sidebar highlights the active solver type based on the
current project.

---

## Implementation Tasks (UI Redesign Phase)

### Step 1 — New `projects` app and data models

- `python manage.py startapp projects`
- Define `Project`, `ProjectMesh`, `ProjectSetup`, `ProjectResult` models
- Add migrations
- Register in `INSTALLED_APPS`

### Step 2 — Base template + sidebar

- Create `templates/base.html` with Bootstrap 5 layout:
  - Fixed left sidebar (250px) listing solver types
  - Top navbar with logo and project name breadcrumb
  - `{% block content %}` main area
- Create `templates/partials/sidebar.html` — solver list with active state

### Step 3 — Home page

- `projects.views.home` — query all `Project` objects, render project cards
- Each card shows: name, solver type, pipeline status (which steps done), last updated
- `[+ New]` button opens the create form

### Step 4 — Create project flow

- `projects.views.create` — form with project name + solver type dropdown
- On submit → create `Project` record → redirect to `/project/<id>/mesh/`

### Step 5 — Pipeline bar partial

- `templates/partials/pipeline_bar.html`
- Three steps: Mesh / Setup / Solve
- Each step is a clickable button if that step or earlier is complete; greyed if locked
- Active step highlighted

### Step 6 — Mesh step

- `projects.views.mesh_step` — mesh type selector + dynamic parameter form
- On submit: generate mesh, save `ProjectMesh`, update `Project.status = 'mesh_done'`
- Preview: generate and show `mesh_plot.png` specific to this project
  - Plot saved as `static/projects/<id>/mesh_plot.png` (per-project paths)
- `[Save & Continue ▶]` → redirect to `/project/<id>/setup/`

### Step 7 — Setup step

- `projects.views.setup_step` — BC form (adapts to solver type and mesh dimension)
- On submit: save `ProjectSetup`, update `Project.status = 'setup_done'`
- `[Save & Continue ▶]` → redirect to `/project/<id>/solve/`

### Step 8 — Solve step

- `projects.views.solve_step` — GET shows results if available; POST runs solver
- On POST: run appropriate solver, save `ProjectResult`
  - Plots saved to `static/projects/<id>/solution_plot.png`
- Update `Project.status = 'solved'`
- Results page shows: mesh plot + solution plot side-by-side, data table, download buttons

### Step 9 — Wire up root URL and redirect old URLs

- Set `path('', include('projects.urls'))` in `FenicsWeb/urls.py`
- Keep old `/SteadyStateThermal/` and `/mesh/` URLs working (for now)

### Step 10 — Per-project static file paths

- Change plot output paths from `static/solution_plot.png` (shared, overwritten)
  to `static/projects/<project_id>/solution_plot.png` (per-project, persistent)
- Update `generate_solution_plot`, `generate_mesh_plot` to accept an output path

---

## Key Design Decisions

### Why a `projects` app instead of modifying existing apps?

The existing `SteadyStateThermal` and `mesh` apps contain solver logic (FEniCS code).
They should stay as **solver libraries** — pure Python functions that take parameters
and return results. The `projects` app is the **orchestration layer** — it owns the
UI, models, and URLs, and calls into the solver libraries.

This separation means adding a new solver in Phase 2 is just:
1. Add `ConvectionDiffusion` solver library (no views/templates)
2. Add `'convection_diffusion'` to `Project.SOLVER_CHOICES`
3. Add a setup form variant for the new parameters
4. The rest of the UI (sidebar, pipeline, project cards) works automatically

### Per-project plot paths

Currently all runs overwrite `static/solution_plot.png`. When a user has multiple
projects, switching between them would show the wrong plots. Fix: save plots to
`static/projects/<id>/solution_plot.png` so each project keeps its own results.

### Session state → database state

The current app passes mesh configuration via Django session (`request.session['mesh_str']`).
This breaks if the session expires. All state moves to the database via `ProjectMesh`,
`ProjectSetup`, and `ProjectResult` — persistent across sessions and browser restarts.

### Sidebar solver availability

Solvers not yet implemented (Phase 2, 3, 6) appear in the sidebar as greyed-out
with a `(coming soon)` label. Clicking them shows a "Not yet available" message.
This makes the roadmap visible to the user without blocking current functionality.

---

## File List to Create / Modify

| File | Action |
| ---- | ------ |
| `projects/__init__.py` | New app |
| `projects/apps.py` | App config |
| `projects/models.py` | Project, ProjectMesh, ProjectSetup, ProjectResult |
| `projects/views.py` | home, create, overview, mesh_step, setup_step, solve_step, delete |
| `projects/urls.py` | URL routing |
| `projects/forms.py` | ProjectForm, MeshStepForm, SetupStepForm |
| `projects/migrations/` | Auto-generated |
| `templates/base.html` | Shared layout with sidebar |
| `templates/partials/sidebar.html` | Solver list |
| `templates/partials/pipeline_bar.html` | Mesh → Setup → Solve indicator |
| `templates/partials/project_card.html` | Reusable project card |
| `templates/projects/home.html` | Main page |
| `templates/projects/create.html` | New project form |
| `templates/projects/overview.html` | Project pipeline view |
| `templates/projects/mesh.html` | Mesh step |
| `templates/projects/setup.html` | Setup step |
| `templates/projects/solve.html` | Results / solve page |
| `SteadyStateThermal/conduction.py` | Accept output_dir for plot paths |
| `mesh/mesh.py` | Accept output_dir for plot paths |
| `FenicsWeb/urls.py` | Add projects URLs as root |
| `FenicsWeb/settings.py` | Add projects to INSTALLED_APPS, set TEMPLATES DIRS |
