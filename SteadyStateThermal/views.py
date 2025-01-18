from django.shortcuts import render, HttpResponse
from .forms import HeatSolverForm, MeshForm
from SteadyStateThermal.conduction import solve_conduction, generate_solution_plot
from mesh.mesh import show_mesh, generate_mesh_plot

# Create your views here.
def conduction(request):
    if request.method == "POST":
        form = HeatSolverForm(request.POST)
        if form.is_valid():
            dimension = form.cleaned_data['dimension']
            left_bc = form.cleaned_data.get('left_bc', 0)
            right_bc = form.cleaned_data.get('right_bc', 0)
            top_bc = form.cleaned_data.get('top_bc', 0)
            bottom_bc = form.cleaned_data.get('bottom_bc', 0)
            front_bc = form.cleaned_data.get('front_bc', 0)
            back_bc = form.cleaned_data.get('back_bc', 0)
            f = form.cleaned_data.get('f', 0)
            k = form.cleaned_data.get('k', 1)

            mesh_str = request.session['mesh_str']
            mesh = show_mesh(mesh_str)
            mesh_plot = generate_mesh_plot(mesh)
            
            solution, solution_list = solve_conduction(dimension, mesh, left_bc, right_bc, top_bc, bottom_bc, front_bc, back_bc, f, k)
            solution_plot = generate_solution_plot(solution)

            context = {
                'mesh_plot': mesh_plot,
                'solution_plot': solution_plot,
                'message': 'Problem Solved',
                'solution_list': solution_list,
            }
            return render(request, 'conduction.html', context)
        else:
            form = MeshForm(request.POST)
            if form.is_valid():
                mesh = form.cleaned_data['mesh']
                request.session['mesh_str'] = mesh.mesh
                mesh = show_mesh(mesh.mesh)
                mesh_plot = generate_mesh_plot(mesh)

                context = {
                    'mesh_plot': mesh_plot,
                    'form': HeatSolverForm(),
                    'message': 'Enter the Boundary Conditions',
                }
                return render(request, 'conduction.html', context)
    else:
        form = MeshForm()
        context = {
            'form': form,
            'message': 'Select the Mesh to Solve',
        }
        return render(request, 'conduction.html', context)