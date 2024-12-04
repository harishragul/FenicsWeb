from django.shortcuts import render, redirect
from django.http import JsonResponse
import mesh.forms as meshform
from mesh.mesh import generate_mesh, generate_mesh_plot, generate_interactive_mesh_plot
    
def generate_mesh_view(request):
    if request.method == "POST":
        try:
            mesh_type = request.session['mesh_type']
            mesh = generate_mesh(mesh_type, request)
            mesh_data = dict(request.POST)
            del mesh_data['csrfmiddlewaretoken']
            mesh_plot = generate_mesh_plot(mesh)
            interactive_mesh_plot = generate_interactive_mesh_plot(mesh)
            context = {
                'mesh_plot': mesh_plot,
                'interactive_mesh_plot': interactive_mesh_plot,
                'mesh_type': mesh_type,
                'message': 'Mesh Generation Completed',
                'data': mesh_data,
            }
            return render(request, 'mesh.html', context=context)
        except Exception as error:
            print(f'-------------------------------{error}----------------------------------')
            mesh_type_form = meshform.MeshType(request.POST)
            if mesh_type_form.is_valid():
                request.session['mesh_type'] = mesh_type_form.cleaned_data['mesh_type']
                mesh_type = request.session['mesh_type']

                if mesh_type == 'interval':
                    mesh_form= meshform.IntervalMesh
                elif mesh_type == 'rectangle':
                    mesh_form= meshform.RectangleMesh
                elif mesh_type == 'box':
                    mesh_form= meshform.BoxMesh
                elif mesh_type == 'unit_interval':
                    mesh_form= meshform.UnitIntervalMesh
                elif mesh_type == 'unit_square':
                    mesh_form= meshform.UnitSquareMesh
                elif mesh_type == 'unit_cube':
                    mesh_form= meshform.UnitCubeMesh
                elif mesh_type == 'circle':
                    mesh_form= meshform.CircleMesh
                elif mesh_type == 'sphere':
                    mesh_form= meshform.SphereMesh

                context = {
                    'form': mesh_form,
                    'mesh_type': mesh_type,
                    'message': 'Enter the Mesh Parameter',
                }
                return render(request, 'mesh.html', context=context)
            

    context = {
        'form': meshform.MeshType,
        'message': 'Select Mesh Type',
    }
    return render(request, 'mesh.html', context=context)

