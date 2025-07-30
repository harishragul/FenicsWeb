from fenics import *
import mesh.forms as meshform

def generate_mesh(mesh_type, request):
    """Generate FEniCS mesh based on the type and parameters."""
    
    def validate_positive_int(value, name):
        """Validate that a value is a positive integer."""
        if value <= 0:
            raise ValueError(f"{name} must be a positive integer, got {value}")
    
    def validate_coordinate_range(x0, x1, axis):
        """Validate that coordinates form a valid range."""
        if x1 <= x0:
            raise ValueError(f"{axis} coordinates invalid: x1 ({x1}) must be greater than x0 ({x0})")
    
    if mesh_type == 'interval':
        data = meshform.IntervalMesh(request.POST)
        if data.is_valid():
            interval_n = data.cleaned_data['interval_n']
            interval_x0 = float(data.cleaned_data['interval_x0'])
            interval_x1 = float(data.cleaned_data['interval_x1'])
            
            # Validation
            validate_positive_int(interval_n, "Number of intervals")
            validate_coordinate_range(interval_x0, interval_x1, "X")
        else:
            raise ValueError("Invalid form data for interval mesh")
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
            
            # Validation
            validate_positive_int(rectangle_nx, "Number of divisions in X")
            validate_positive_int(rectangle_ny, "Number of divisions in Y")
            validate_coordinate_range(rectangle_x0, rectangle_x1, "X")
            validate_coordinate_range(rectangle_y0, rectangle_y1, "Y")
        else:
            raise ValueError("Invalid form data for rectangle mesh")
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

def generate_mesh_plot(mesh):
    plot(mesh)

    if not os.path.exists(STATIC_DIR):
        os.makedirs(STATIC_DIR)

    plot_filename = os.path.join(STATIC_DIR, 'mesh_plot.png')
    plt.savefig(plot_filename)
    plt.close()

    return 'mesh_plot.png'

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