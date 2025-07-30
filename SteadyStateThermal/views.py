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
            f = form.cleaned_data.get('f', 0)
            k = form.cleaned_data.get('k', 1)
            
            # Build boundary condition data structure
            bc_data = {}
            
            # Define boundary names based on dimension
            if dimension == '1D':
                boundary_names = ['left', 'right']
            elif dimension == '2D':
                boundary_names = ['left', 'right', 'top', 'bottom']
            elif dimension == '3D':
                boundary_names = ['left', 'right', 'top', 'bottom', 'front', 'back']
            
            # Process each boundary
            for boundary in boundary_names:
                bc_type = form.cleaned_data.get(f'{boundary}_bc_type')
                if bc_type:
                    bc_info = {'type': bc_type}
                    
                    if bc_type == 'dirichlet':
                        bc_value = form.cleaned_data.get(f'{boundary}_bc')
                        if bc_value is not None:
                            bc_info['value'] = bc_value
                        else:
                            # Skip this boundary if no value provided for Dirichlet BC
                            continue
                    elif bc_type == 'convection':
                        h = form.cleaned_data.get(f'{boundary}_h')
                        t_amb = form.cleaned_data.get(f'{boundary}_t_amb')
                        if h is not None and t_amb is not None:
                            bc_info['h'] = h
                            bc_info['t_amb'] = t_amb
                        else:
                            # Skip this boundary if convection parameters are missing
                            continue
                    
                    bc_data[boundary] = bc_info

            mesh_str = request.session['mesh_str']
            mesh = show_mesh(mesh_str)
            mesh_plot = generate_mesh_plot(mesh)
            
            try:
                solution, solution_list = solve_conduction(dimension, mesh, bc_data, f, k)
                solution_plots = generate_solution_plot(solution)

                context = {
                    'mesh_plot': mesh_plot,
                    'solution_plot_interactive': solution_plots['interactive'],
                    'solution_plot_static': solution_plots['static'],
                    'message': 'Problem Solved Successfully',
                    'solution_list': solution_list,
                }
                return render(request, 'conduction.html', context)
            except Exception as e:
                context = {
                    'mesh_plot': mesh_plot,
                    'message': f'Error solving problem: {str(e)}',
                    'form': HeatSolverForm(request.POST),
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