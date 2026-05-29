from django.shortcuts import render
from .forms import HeatSolverForm, MeshForm
from SteadyStateThermal.conduction import solve_conduction, generate_solution_plot
from mesh.mesh import show_mesh, generate_mesh_plot


def conduction(request):
    if request.method == "POST":
        form = HeatSolverForm(request.POST)
        if form.is_valid():
            d = form.cleaned_data
            dimension = d['dimension']

            bc_types = {
                'left':   d['left_bc_type'],
                'right':  d['right_bc_type'],
                'bottom': d['bottom_bc_type'],
                'top':    d['top_bc_type'],
                'front':  d['front_bc_type'],
                'back':   d['back_bc_type'],
            }
            neumann_fluxes = {
                'left':   d.get('left_bc') or 0.0,
                'right':  d.get('right_bc') or 0.0,
                'bottom': d.get('bottom_bc') or 0.0,
                'top':    d.get('top_bc') or 0.0,
                'front':  d.get('front_bc') or 0.0,
                'back':   d.get('back_bc') or 0.0,
            }
            robin_hs = {
                'left':   d.get('left_h') or 0.0,
                'right':  d.get('right_h') or 0.0,
                'bottom': d.get('bottom_h') or 0.0,
                'top':    d.get('top_h') or 0.0,
                'front':  d.get('front_h') or 0.0,
                'back':   d.get('back_h') or 0.0,
            }
            robin_u_infs = {
                'left':   d.get('left_u_inf') or 0.0,
                'right':  d.get('right_u_inf') or 0.0,
                'bottom': d.get('bottom_u_inf') or 0.0,
                'top':    d.get('top_u_inf') or 0.0,
                'front':  d.get('front_u_inf') or 0.0,
                'back':   d.get('back_u_inf') or 0.0,
            }

            mesh_str = request.session['mesh_str']
            mesh = show_mesh(mesh_str)
            mesh_plot = generate_mesh_plot(mesh)

            u_sol, solution_list = solve_conduction(
                dimension, mesh,
                left_bc=d.get('left_bc') or 0.0,
                right_bc=d.get('right_bc') or 0.0,
                top_bc=d.get('top_bc') or 0.0,
                bottom_bc=d.get('bottom_bc') or 0.0,
                front_bc=d.get('front_bc') or 0.0,
                back_bc=d.get('back_bc') or 0.0,
                f=d.get('f', 0.0),
                k=d.get('k', 1.0),
                bc_types=bc_types,
                neumann_fluxes=neumann_fluxes,
                robin_hs=robin_hs,
                robin_u_infs=robin_u_infs,
            )
            solution_plot = generate_solution_plot(u_sol, mesh)

            context = {
                'mesh_plot': mesh_plot,
                'solution_plot': solution_plot,
                'message': 'Problem Solved',
                'solution_list': solution_list,
            }
            return render(request, 'conduction.html', context)

        else:
            mesh_form = MeshForm(request.POST)
            if mesh_form.is_valid():
                mesh_obj = mesh_form.cleaned_data['mesh']
                request.session['mesh_str'] = mesh_obj.mesh
                mesh = show_mesh(mesh_obj.mesh)
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


def validation(request):
    """Run h-refinement convergence study for 1D steady-state conduction."""
    from fenics import UnitIntervalMesh, errornorm, Expression
    import math

    results = []
    ns = [4, 8, 16, 32, 64]
    u_exact_expr = Expression("x[0]*(1 - x[0])", degree=4)

    prev_err = None
    for n in ns:
        mesh = UnitIntervalMesh(n)
        u_sol, _ = solve_conduction(
            '1D', mesh,
            left_bc=0, right_bc=0,
            top_bc=0, bottom_bc=0, front_bc=0, back_bc=0,
            f=2, k=1,
        )
        h = 1.0 / n
        err = errornorm(u_exact_expr, u_sol, 'L2')
        rate = math.log(prev_err / err) / math.log(2) if prev_err else '—'
        results.append({
            'n': n,
            'h': f'{h:.4f}',
            'l2_error': f'{err:.2e}',
            'rate': f'{rate:.2f}' if isinstance(rate, float) else rate,
        })
        prev_err = err

    context = {
        'results': results,
        'exact': 'u(x) = x(1 − x)',
        'problem': '−u″ = 2,  u(0) = u(1) = 0',
        'expected_rate': '2.0 (P1 elements, L2 norm)',
    }
    return render(request, 'validation.html', context)
