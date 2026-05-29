import os
from django.shortcuts import render, redirect, get_object_or_404
from django.conf import settings

from .models import Project, ProjectMesh, ProjectSetup, ProjectResult
from .forms import ProjectCreateForm, MeshStepForm, SetupStepForm
from .solver_bridge import run_project_solver
from mesh.mesh import build_mesh_from_params, generate_mesh_plot
from SteadyStateThermal.conduction import generate_solution_plot

FACES_1D = ['left', 'right']
FACES_2D = ['left', 'right', 'bottom', 'top']
FACES_3D = ['left', 'right', 'bottom', 'top', 'front', 'back']

MESH_TYPE_TO_DIM = {
    'interval':      '1D',
    'unit_interval': '1D',
    'rectangle':     '2D',
    'unit_square':   '2D',
    'box':           '3D',
    'unit_cube':     '3D',
}


def _dimension_from_mesh(project):
    """Return '1D'/'2D'/'3D' derived from the saved mesh type."""
    try:
        mt = project.mesh_config.mesh_params.get('mesh_type', '')
        return MESH_TYPE_TO_DIM.get(mt, '3D')
    except Exception:
        return '3D'


# ── Helpers ───────────────────────────────────────────────────────────────────

def _project_media_dir(project):
    path = os.path.join(settings.MEDIA_ROOT, 'projects', str(project.pk))
    os.makedirs(path, exist_ok=True)
    return path


def _media_rel(project, filename):
    return os.path.join('projects', str(project.pk), filename)


def _build_mesh_from_form(form):
    """Return (fenics_mesh, mesh_str, mesh_params_dict) from validated MeshStepForm."""
    d  = form.cleaned_data
    mt = d['mesh_type']
    if mt == 'interval':
        params = {'n': d['n'], 'x0': d['x0'], 'x1': d['x1']}
    elif mt == 'unit_interval':
        params = {'n': d['n']}
    elif mt == 'rectangle':
        params = {'x0': d['x0'], 'y0': d['y0'], 'x1': d['x1'], 'y1': d['y1'],
                  'nx': d['nx'], 'ny': d['ny']}
    elif mt == 'unit_square':
        params = {'nx': d['nx'], 'ny': d['ny']}
    elif mt == 'box':
        params = {'x0': d['x0'], 'y0': d['y0'], 'z0': d['z0'],
                  'x1': d['x1'], 'y1': d['y1'], 'z1': d['z1'],
                  'nx': d['nx'], 'ny': d['ny'], 'nz': d['nz']}
    elif mt == 'unit_cube':
        params = {'nx': d['nx'], 'ny': d['ny'], 'nz': d['nz']}
    else:
        raise ValueError(f'Unknown mesh type: {mt}')
    mesh, mesh_str = build_mesh_from_params(mt, params)
    return mesh, mesh_str, {'mesh_type': mt, **params}


def _extract_setup(form, dimension):
    """Convert SetupStepForm → (bc_config dict, solver_params dict)."""
    d  = form.cleaned_data
    bc = {}
    for face in FACES_3D:
        ftype = d.get(f'{face}_type', 'dirichlet')
        entry = {'type': ftype}
        if ftype in ('dirichlet', 'neumann'):
            entry['value'] = float(d.get(f'{face}_value') or 0.0)
        elif ftype == 'robin':
            entry['h']     = float(d.get(f'{face}_h')     or 0.0)
            entry['u_inf'] = float(d.get(f'{face}_u_inf') or 0.0)
        bc[face] = entry
    solver_params = {'k': float(d['k']), 'f': float(d['f'])}
    return bc, solver_params


def _setup_to_initial(setup):
    """Convert saved ProjectSetup → form initial dict."""
    initial = {
        'k': setup.solver_params.get('k', 1.0),
        'f': setup.solver_params.get('f', 0.0),
    }
    for face, cfg in setup.bc_config.items():
        initial[f'{face}_type']  = cfg.get('type', 'dirichlet')
        initial[f'{face}_value'] = cfg.get('value', 0.0)
        initial[f'{face}_h']     = cfg.get('h', 0.0)
        initial[f'{face}_u_inf'] = cfg.get('u_inf', 0.0)
    return initial


def _active_faces(dimension):
    if dimension == '1D':
        return FACES_1D
    elif dimension == '2D':
        return FACES_2D
    return FACES_3D


# ── Views ─────────────────────────────────────────────────────────────────────

def home(request):
    projects = Project.objects.all().order_by('-updated_at')
    return render(request, 'projects/home.html', {'projects': projects})


def create(request):
    solver_preset = request.GET.get('solver', '')
    initial = {'solver_type': solver_preset} if solver_preset else {}
    if request.method == 'POST':
        form = ProjectCreateForm(request.POST)
        if form.is_valid():
            project = form.save()
            return redirect('project_mesh', pk=project.pk)
    else:
        form = ProjectCreateForm(initial=initial)
    return render(request, 'projects/create.html', {'form': form})


def overview(request, pk):
    project = get_object_or_404(Project, pk=pk)
    return render(request, 'projects/overview.html', {'project': project})


def mesh_step(request, pk):
    project = get_object_or_404(Project, pk=pk)
    if request.method == 'POST':
        form = MeshStepForm(request.POST)
        if form.is_valid():
            try:
                mesh, mesh_str, mesh_params = _build_mesh_from_form(form)
                output_dir   = _project_media_dir(project)
                generate_mesh_plot(mesh, output_dir=output_dir)
                plot_rel     = _media_rel(project, 'mesh_plot.png')
                ProjectMesh.objects.update_or_create(
                    project=project,
                    defaults={
                        'mesh_type':      form.cleaned_data['mesh_type'],
                        'mesh_str':       mesh_str,
                        'mesh_params':    mesh_params,
                        'mesh_plot_path': plot_rel,
                    },
                )
                if project.status == 'new':
                    project.status = 'mesh_done'
                    project.save()
                return redirect('project_setup', pk=project.pk)
            except Exception as exc:
                form.add_error(None, str(exc))
    else:
        initial = {}
        if hasattr(project, 'mesh_config'):
            p = project.mesh_config.mesh_params
            initial = {k: v for k, v in p.items() if k != 'mesh_type'}
            initial['mesh_type'] = p.get('mesh_type', 'interval')
        form = MeshStepForm(initial=initial)
    return render(request, 'projects/mesh.html', {'project': project, 'form': form})


def setup_step(request, pk):
    project = get_object_or_404(Project, pk=pk)
    # Dimension is fixed by the mesh — not a user choice.
    dimension = _dimension_from_mesh(project)
    active_faces = {'1D': FACES_1D, '2D': FACES_2D, '3D': FACES_3D}[dimension]

    if request.method == 'POST':
        form = SetupStepForm(request.POST)
        if form.is_valid():
            bc_config, solver_params = _extract_setup(form, dimension)
            ProjectSetup.objects.update_or_create(
                project=project,
                defaults={
                    'dimension':     dimension,
                    'bc_config':     bc_config,
                    'solver_params': solver_params,
                },
            )
            if project.status in ('new', 'mesh_done'):
                project.status = 'setup_done'
            project.save()
            return redirect('project_overview', pk=project.pk)
    else:
        initial = {}
        if hasattr(project, 'setup_config'):
            initial = _setup_to_initial(project.setup_config)
        form = SetupStepForm(initial=initial)

    face_groups = [
        {
            'name':        face,
            'label':       f'{face.capitalize()} face',
            'type_field':  form[f'{face}_type'],
            'value_field': form[f'{face}_value'],
            'h_field':     form[f'{face}_h'],
            'u_inf_field': form[f'{face}_u_inf'],
        }
        for face in active_faces   # only the faces that match this mesh's dimension
    ]
    return render(request, 'projects/setup.html', {
        'project':     project,
        'form':        form,
        'face_groups': face_groups,
        'dimension':   dimension,
    })


def solve_step(request, pk):
    project = get_object_or_404(Project, pk=pk)
    if request.method == 'POST':
        try:
            u_sol, solution_list, elapsed = run_project_solver(project)
            from mesh.mesh import show_mesh
            mesh       = show_mesh(project.mesh_config.mesh_str)
            output_dir = _project_media_dir(project)
            generate_solution_plot(u_sol, mesh, output_dir=output_dir)
            generate_mesh_plot(mesh, output_dir=output_dir)
            ProjectResult.objects.update_or_create(
                project=project,
                defaults={
                    'mesh_plot_path':     _media_rel(project, 'mesh_plot.png'),
                    'solution_plot_path': _media_rel(project, 'solution_plot.png'),
                    'solution_data':      solution_list,
                    'solve_time_s':       elapsed,
                },
            )
            project.status = 'solved'
            project.save()
        except Exception as exc:
            return render(request, 'projects/solve.html', {
                'project': project,
                'error':   str(exc),
            })
        return redirect('project_solve', pk=project.pk)
    return render(request, 'projects/solve.html', {'project': project})


def delete(request, pk):
    project = get_object_or_404(Project, pk=pk)
    if request.method == 'POST':
        project.delete()
        return redirect('home')
    return redirect('project_overview', pk=pk)
