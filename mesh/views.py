from django.shortcuts import render, HttpResponse
from .forms import *
from fenics import *

# Create your views here.
def generate_mesh(request):
    if request.POST.get('mesh_type'):
        mesh_type = request.POST.get('mesh_type')

        context = {
            'type': mesh_type,
        }

        if mesh_type == 'interval':
            context['form'] = IntervalMesh
        elif mesh_type == 'rectangle':
            context['form'] = RectangleMesh
        elif mesh_type == 'box':
            context['form'] = BoxMesh
        elif mesh_type == 'unit_interval':
            context['form'] = UnitIntervalMesh
        elif mesh_type == 'unit_square':
            context['form'] = UnitSquareMesh
        elif mesh_type == 'unit_cube':
            context['form'] = UnitCubeMesh
        elif mesh_type == 'circle':
            context['form'] = CircleMesh
        elif mesh_type == 'sphere':
            context['form'] = SphereMesh
        elif mesh_type == 'custom':
            context['form'] = CustomMesh

        return render(request, 'generate_mesh.html', context)

    else:
        context = {
            'form': MeshType,
        }

        return render(request, 'generate_mesh.html', context)
    
def show_mesh(request):
    if request.POST.get('mesh_type'):
        mesh_type = request.POST.get('mesh_type')
        mesh_data = request.POST.items()

        if mesh_type == 'interval':
            mesh = IntervalMesh(mesh_data['interval_n'], mesh_data['interval_x0'], mesh_data['interval_x1'])
        elif mesh_type == 'rectangle':
            mesh = RectangleMesh
        elif mesh_type == 'box':
            mesh = BoxMesh
        elif mesh_type == 'unit_interval':
            mesh = UnitIntervalMesh
        elif mesh_type == 'unit_square':
            mesh = UnitSquareMesh
        elif mesh_type == 'unit_cube':
            mesh = UnitCubeMesh
        elif mesh_type == 'circle':
            mesh = CircleMesh
        elif mesh_type == 'sphere':
            mesh = SphereMesh
        elif mesh_type == 'custom':
            mesh = CustomMesh
