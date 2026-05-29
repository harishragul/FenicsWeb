from fenics import *
import mesh.forms as meshform

def generate_mesh(mesh_type, request):
    """Generate FEniCS mesh based on the type and parameters."""
    if mesh_type == 'interval':
        data = meshform.IntervalMesh(request.POST)
        if data.is_valid():
            interval_n = data.cleaned_data['interval_n']
            interval_x0 = float(data.cleaned_data['interval_x0'])
            interval_x1 = float(data.cleaned_data['interval_x1'])
        else:
            interval_n = 0
            interval_x0 = 0.0
            interval_x1 = 0.0
        mesh_str = f"{mesh_type},{interval_n},{interval_x0},{interval_x1}"
        return IntervalMesh(interval_n, interval_x0, interval_x1), mesh_str
    
    elif mesh_type == 'rectangle':
        data = meshform.RectangleMesh(request.POST)
        if data.is_valid():
            rectangle_x0 = data.cleaned_data['rectangle_x0']
            rectangle_y0 = data.cleaned_data['rectangle_y0']
            rectangle_x1 = data.cleaned_data['rectangle_x1']
            rectangle_y1 = data.cleaned_data['rectangle_y1']
            rectangle_nx = data.cleaned_data['rectangle_nx']
            rectangle_ny = data.cleaned_data['rectangle_ny']
        else:
            rectangle_x0 = 0
            rectangle_y0 = 0
            rectangle_x1 = 0
            rectangle_y1 = 0
            rectangle_nx = 0
            rectangle_ny = 0
        mesh_str = f"{mesh_type},{rectangle_x0},{rectangle_y0},{rectangle_x1},{rectangle_y1},{rectangle_nx},{rectangle_ny}"
        return RectangleMesh(Point(rectangle_x0, rectangle_y0), Point(rectangle_x1, rectangle_y1), rectangle_nx, rectangle_ny), mesh_str
    
    elif mesh_type == 'box':
        data = meshform.BoxMesh(request.POST)
        if data.is_valid():
            box_x0 = data.cleaned_data['box_x0']
            box_y0 = data.cleaned_data['box_y0']
            box_z0 = data.cleaned_data['box_z0']
            box_x1 = data.cleaned_data['box_x1']
            box_y1 = data.cleaned_data['box_y1']
            box_z1 = data.cleaned_data['box_z1']
            box_nx = data.cleaned_data['box_nx']
            box_ny = data.cleaned_data['box_ny']
            box_nz = data.cleaned_data['box_nz']
        else:
            box_x0 = 0
            box_y0 = 0
            box_z0 = 0
            box_x1 = 0
            box_y1 = 0
            box_z1 = 0
            box_nx = 0
            box_ny = 0
            box_nz = 0
        mesh_str = f"{mesh_type},{box_x0},{box_y0},{box_z0},{box_x1},{box_y1},{box_z1},{box_nx},{box_ny},{box_nz}"
        return BoxMesh(Point(box_x0, box_y0, box_z0), Point(box_x1, box_y1, box_z1), box_nx, box_ny, box_nz), mesh_str
    
    elif mesh_type == 'unit_interval':
        data = meshform.UnitIntervalMesh(request.POST)
        if data.is_valid():
            unit_nx = data.cleaned_data['unit_nx']
        else:
            unit_nx = 0
        mesh_str = f'{mesh_type},{unit_nx}'
        return UnitIntervalMesh(unit_nx), mesh_str
    
    elif mesh_type == 'unit_square':
        data = meshform.UnitSquareMesh(request.POST)
        if data.is_valid():
            unit_nx = data.cleaned_data['unit_nx']
            unit_ny = data.cleaned_data['unit_ny']
        else:
            unit_nx = 0
            unit_ny = 0
        mesh_str = f'{mesh_type},{unit_nx},{unit_ny}'
        return UnitSquareMesh(unit_nx, unit_ny), mesh_str
    
    elif mesh_type == 'unit_cube':
        data = meshform.UnitCubeMesh(request.POST)

        if data.is_valid():
            unit_nx = data.cleaned_data['unit_nx']
            unit_ny = data.cleaned_data['unit_ny']
            unit_nz = data.cleaned_data['unit_nz']
        else:
            unit_nx = 0
            unit_ny = 0
            unit_nz = 0
        mesh_str = f'{mesh_type},{unit_nx},{unit_ny},{unit_nz}'
        return UnitCubeMesh(unit_nx, unit_ny, unit_nz), mesh_str
    
    else:
        raise ValueError("Invalid mesh type.")
    
#Show Mesh Function
def show_mesh(mesh_str):
    mesh_list = mesh_str.split(',')

    if mesh_list[0] == 'interval':
        interval_n = int(mesh_list[1])
        interval_x0 = float(mesh_list[2])
        interval_x1 = float(mesh_list[3])
        return IntervalMesh(interval_n, interval_x0, interval_x1)
    
    elif mesh_list[0] == 'rectangle':
        rectangle_x0 = float(mesh_list[1])
        rectangle_y0 = float(mesh_list[2])
        rectangle_x1 = float(mesh_list[3])
        rectangle_y1 = float(mesh_list[4])
        rectangle_nx = int(mesh_list[5])
        rectangle_ny = int(mesh_list[6])
        return RectangleMesh(Point(rectangle_x0, rectangle_y0), Point(rectangle_x1, rectangle_y1), rectangle_nx, rectangle_ny)
    
    elif mesh_list[0] == 'box':
        box_x0 = float(mesh_list[1])
        box_y0 = float(mesh_list[2])
        box_z0 = float(mesh_list[3])
        box_x1 = float(mesh_list[4])
        box_y1 = float(mesh_list[5])
        box_z1 = float(mesh_list[6])
        box_nx = int(mesh_list[7])
        box_ny = int(mesh_list[8])
        box_nz = int(mesh_list[9])
        return BoxMesh(Point(box_x0, box_y0, box_z0), Point(box_x1, box_y1, box_z1), box_nx, box_ny, box_nz)
    
    elif mesh_list[0] == 'unit_interval':
        unit_nx = int(mesh_list[1])
        return UnitIntervalMesh(unit_nx)
    
    elif mesh_list[0] == 'unit_square':
        unit_nx = int(mesh_list[1])
        unit_ny = int(mesh_list[2])
        return UnitSquareMesh(unit_nx, unit_ny)
    
    elif mesh_list[0] == 'unit_cube':
        unit_nx = int(mesh_list[1])
        unit_ny = int(mesh_list[2])
        unit_nz = int(mesh_list[3])
        return UnitCubeMesh(unit_nx, unit_ny, unit_nz)
    else:
        raise ValueError("Invalid mesh type.")
    

import os
import matplotlib.pyplot as plt
from FenicsWeb.settings import BASE_DIR
import matplotlib
matplotlib.use("Agg")

STATIC_DIR = os.path.join(BASE_DIR, 'static') 

def build_mesh_from_params(mesh_type, params):
    """Create a FEniCS mesh from a plain dict — no Django request needed."""
    if mesh_type == 'interval':
        mesh = IntervalMesh(params['n'], params['x0'], params['x1'])
        mesh_str = f"interval,{params['n']},{params['x0']},{params['x1']}"
    elif mesh_type == 'unit_interval':
        mesh = UnitIntervalMesh(params['n'])
        mesh_str = f"unit_interval,{params['n']}"
    elif mesh_type == 'rectangle':
        mesh = RectangleMesh(
            Point(params['x0'], params['y0']),
            Point(params['x1'], params['y1']),
            params['nx'], params['ny'],
        )
        mesh_str = (f"rectangle,{params['x0']},{params['y0']},"
                    f"{params['x1']},{params['y1']},{params['nx']},{params['ny']}")
    elif mesh_type == 'unit_square':
        mesh = UnitSquareMesh(params['nx'], params['ny'])
        mesh_str = f"unit_square,{params['nx']},{params['ny']}"
    elif mesh_type == 'box':
        mesh = BoxMesh(
            Point(params['x0'], params['y0'], params['z0']),
            Point(params['x1'], params['y1'], params['z1']),
            params['nx'], params['ny'], params['nz'],
        )
        mesh_str = (f"box,{params['x0']},{params['y0']},{params['z0']},"
                    f"{params['x1']},{params['y1']},{params['z1']},"
                    f"{params['nx']},{params['ny']},{params['nz']}")
    elif mesh_type == 'unit_cube':
        mesh = UnitCubeMesh(params['nx'], params['ny'], params['nz'])
        mesh_str = f"unit_cube,{params['nx']},{params['ny']},{params['nz']}"
    else:
        raise ValueError(f"Unknown mesh type: {mesh_type}")
    return mesh, mesh_str


def generate_mesh_plot(mesh, output_dir=None):
    """Plot mesh topology using direct matplotlib (avoids FEniCS gca() compat issue)."""
    import matplotlib.tri as mtri
    import numpy as np

    coords = mesh.coordinates()
    cells = mesh.cells()
    dim = mesh.geometry().dim()

    target_dir = output_dir if output_dir is not None else STATIC_DIR
    os.makedirs(target_dir, exist_ok=True)
    plot_filename = os.path.join(target_dir, 'mesh_plot.png')

    if dim == 1:
        fig, ax = plt.subplots(figsize=(8, 3))
        x = coords.flatten()
        ax.scatter(x, np.zeros_like(x), s=15, color='steelblue', zorder=3)
        for cell in cells:
            ax.plot(coords[cell, 0], [0, 0], 'steelblue', linewidth=1.5)
        ax.set_xlabel('x')
        ax.set_yticks([])
        ax.set_title(f'1D Mesh  ({mesh.num_cells()} cells)')

    elif dim == 2:
        fig, ax = plt.subplots(figsize=(7, 6))
        triang = mtri.Triangulation(coords[:, 0], coords[:, 1], cells)
        ax.triplot(triang, color='steelblue', linewidth=0.6)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_aspect('equal')
        ax.set_title(f'2D Mesh  ({mesh.num_cells()} cells)')

    elif dim == 3:
        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111, projection='3d')
        plotted = 0
        for cell in cells:
            verts = coords[cell]
            n = len(cell)
            for i in range(n):
                for j in range(i + 1, n):
                    ax.plot(*zip(verts[i], verts[j]),
                            color='steelblue', linewidth=0.3, alpha=0.4)
            plotted += 1
            if plotted > 500:   # cap edge drawing for large meshes
                break
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        ax.set_title(f'3D Mesh  ({mesh.num_cells()} cells)')

    plt.tight_layout()
    plt.savefig(plot_filename, dpi=100, bbox_inches='tight')
    plt.close()
    return 'mesh_plot.png'




import plotly.graph_objects as go
import json, plotly
import plotly.utils
import json

def generate_interactive_mesh_plot(mesh, u_sol=None):
    """Return a Plotly Mesh3d figure as JSON. Colours vertices by u_sol if provided."""
    coords = mesh.coordinates()
    connectivity = mesh.cells()

    intensity = u_sol.compute_vertex_values(mesh).tolist() if u_sol is not None else None

    mesh3d_kwargs = dict(
        x=coords[:, 0].tolist(),
        y=coords[:, 1].tolist(),
        z=coords[:, 2].tolist() if mesh.geometry().dim() == 3 else [0] * len(coords),
        i=connectivity[:, 0].tolist(),
        j=connectivity[:, 1].tolist(),
        k=connectivity[:, 2].tolist(),
        opacity=0.8,
    )
    if intensity is not None:
        mesh3d_kwargs['intensity'] = intensity
        mesh3d_kwargs['colorscale'] = 'Viridis'
        mesh3d_kwargs['colorbar'] = {'title': 'Temperature (K)'}

    fig = go.Figure(data=[go.Mesh3d(**mesh3d_kwargs)])
    fig.update_layout(margin=dict(l=0, r=0, b=0, t=30))
    return json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)