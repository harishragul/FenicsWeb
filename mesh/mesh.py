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
        mesh_str = f"{mesh_str},{box_x0},{box_y0},{box_z0},{box_x1},{box_y1},{box_z1},{box_nx},{box_ny},{box_nz}"
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
    
    elif mesh_type == 'circle':
        data = meshform.CircleMesh(request.POST)
        if data.is_valid():
            circle_xc = data.cleaned_data['circle_xc']
            circle_yc = data.cleaned_data['circle_yc']
            circle_radius = data.cleaned_data['circle_radius']
            circle_resolution = data.cleaned_data['circle_resolution']
        else:
            circle_xc = 0
            circle_yc = 0
            circle_radius = 0
            circle_resolution = 0
        mesh_str = f'{mesh_type},{circle_xc},{circle_yc},{circle_radius},{circle_resolution}'
        return CircleMesh(Point(circle_xc, circle_yc),circle_radius,circle_resolution), mesh_str
    
    elif mesh_type == 'sphere':
        data = meshform.SphereMesh(request.POST)
        if data.is_valid():
            sphere_xc = data.cleaned_data['sphere_xc']
            sphere_yc = data.cleaned_data['sphere_yc']
            sphere_zc = data.cleaned_data['sphere_zc']
            sphere_radius = data.cleaned_data['sphere_radius']
            sphere_resolution = data.cleaned_data['sphere_resolution']
        else:
            sphere_xc = 0
            sphere_yc = 0
            sphere_zc = 0
            sphere_radius = 0
            sphere_resolution = 0
            mesh_str = f'{mesh_type},{sphere_xc},{sphere_yc},{sphere_zc},{sphere_radius},{sphere_resolution}'
        return SphereMesh(Point(sphere_xc, sphere_yc, sphere_zc), sphere_radius,sphere_resolution), mesh_str
    else:
        raise ValueError("Invalid mesh type.")
    

import os
import matplotlib.pyplot as plt
from FenicsWeb.settings import BASE_DIR
import matplotlib
matplotlib.use("Agg")

STATIC_DIR = os.path.join(BASE_DIR, 'static') 

def generate_mesh_plot(mesh):
    plot(mesh)

    if not os.path.exists(STATIC_DIR):
        os.makedirs(STATIC_DIR)

    plot_filename = os.path.join(STATIC_DIR, 'temp_plot.png')
    plt.savefig(plot_filename)
    plt.close()

    return 'temp_plot.png'

import plotly.graph_objects as go
import json, plotly
import plotly.utils
import json

def generate_interactive_mesh_plot(mesh):
    coords = mesh.coordinates()
    connectivity = mesh.cells()

    fig = go.Figure(data=[go.Mesh3d(
        x=coords[:, 0], y=coords[:, 1], z=coords[:, 2] if mesh.geometry().dim() == 3 else [0] * len(coords),
        i=connectivity[:, 0], j=connectivity[:, 1], k=connectivity[:, 2]
    )])

    return json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)