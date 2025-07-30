from django.shortcuts import render, redirect
from django.http import JsonResponse
import mesh.forms as meshform
from mesh.mesh import generate_mesh, generate_mesh_plot, generate_interactive_mesh_plot
from mesh.models import Mesh

def generate_mesh_view(request):
    if request.method == "POST":
        try:
            # Check if mesh_type is in session
            if 'mesh_type' not in request.session:
                # Try to get mesh type from current form
                mesh_type_form = meshform.MeshType(request.POST)
                if mesh_type_form.is_valid():
                    request.session['mesh_type'] = mesh_type_form.cleaned_data['mesh_type']
                    request.session['mesh_name'] = mesh_type_form.cleaned_data['mesh_name']
                else:
                    context = {
                        'form': meshform.MeshType,
                        'message': 'Please select a valid mesh type',
                        'error': True
                    }
                    return render(request, 'mesh.html', context=context)
            
            mesh_type = request.session['mesh_type']
            mesh, mesh_str = generate_mesh(mesh_type, request)
            Mesh(name=request.session['mesh_name'], mesh=mesh_str).save()
            
            mesh_data = dict(request.POST)
            del mesh_data['csrfmiddlewaretoken']
            mesh_plot = generate_mesh_plot(mesh)
            #interactive_mesh_plot = generate_interactive_mesh_plot(mesh)
            context = {
                'mesh_plot': mesh_plot,
                #'interactive_mesh_plot': interactive_mesh_plot,
                'mesh_type': mesh_type,
                'message': 'Mesh Generation Completed',
                'data': mesh_data,
            }
            return render(request, 'mesh.html', context=context)
            
        except ValueError as e:
            # Handle validation errors specifically
            context = {
                'form': meshform.MeshType,
                'message': f'Input validation error: {str(e)}',
                'error': True
            }
            return render(request, 'mesh.html', context=context)
            
        except Exception as error:
            # Log the actual error for debugging (in production, use proper logging)
            print(f"Mesh generation error: {str(error)}")
            
            # Try to handle mesh type form submission
            mesh_type_form = meshform.MeshType(request.POST)
            if mesh_type_form.is_valid():
                request.session['mesh_type'] = mesh_type_form.cleaned_data['mesh_type']
                request.session['mesh_name'] = mesh_type_form.cleaned_data['mesh_name']
                mesh_type = request.session['mesh_type']

                # Get appropriate form for mesh parameters
                mesh_form_mapping = {
                    'interval': meshform.IntervalMesh,
                    'rectangle': meshform.RectangleMesh,
                    'box': meshform.BoxMesh,
                    'unit_interval': meshform.UnitIntervalMesh,
                    'unit_square': meshform.UnitSquareMesh,
                    'unit_cube': meshform.UnitCubeMesh
                }
                
                mesh_form = mesh_form_mapping.get(mesh_type)
                if mesh_form:
                    context = {
                        'form': mesh_form,
                        'mesh_type': mesh_type,
                        'message': f'Enter the Mesh Parameters for {mesh_type}',
                    }
                    return render(request, 'mesh.html', context=context)
            
            # Fallback error case
            context = {
                'form': meshform.MeshType,
                'message': 'An error occurred during mesh generation. Please try again.',
                'error': True
            }
            return render(request, 'mesh.html', context=context)
            

    context = {
        'form': meshform.MeshType,
        'message': 'Select Mesh Type',
    }
    return render(request, 'mesh.html', context=context)

