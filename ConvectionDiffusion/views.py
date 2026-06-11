from django.shortcuts import render
from fenics import UnitIntervalMesh, UnitSquareMesh

from .forms import ConvectionDiffusionForm
from .solver import (
    solve_convection_diffusion, compute_peclet_number,
    generate_convection_solution_plot, peclet_study_1d, generate_peclet_plot,
)

PE_STUDY_VALUES = [1, 10, 50, 100]


def convection(request):
    if request.method == "POST":
        form = ConvectionDiffusionForm(request.POST)
        if form.is_valid():
            d = form.cleaned_data
            dimension = d['dimension']
            n = d['n']

            mesh = UnitIntervalMesh(n) if dimension == '1D' else UnitSquareMesh(n, n)

            methods = ['galerkin', 'supg'] if d['method'] == 'both' else [d['method']]
            results = []
            for method in methods:
                u_sol, solution_list = solve_convection_diffusion(
                    dimension, mesh,
                    k=d['k'], bx=d['bx'], by=d['by'], f=d['f'],
                    left_bc=d['left_bc'], right_bc=d['right_bc'],
                    bottom_bc=d['bottom_bc'], top_bc=d['top_bc'],
                    method=method,
                )
                plot = generate_convection_solution_plot(
                    u_sol, mesh, d['bx'], d['by'], label=method,
                )
                results.append({
                    'method': method,
                    'plot': plot,
                    'solution_list': solution_list,
                })

            peclet = compute_peclet_number(mesh, d['k'], d['bx'], d['by'])

            context = {
                'form': form,
                'results': results,
                'peclet': peclet,
                'message': 'Problem Solved',
            }
            return render(request, 'convection.html', context)
    else:
        form = ConvectionDiffusionForm()

    return render(request, 'convection.html', {
        'form': form,
        'message': 'Enter Convection-Diffusion Parameters',
    })


def peclet_study(request):
    """1D Peclet number study comparing standard Galerkin and SUPG against
    the exact solution of -u'' + Pe*u' = 0, u(0)=0, u(1)=1.
    """
    results = peclet_study_1d(PE_STUDY_VALUES)
    plot = generate_peclet_plot(results)
    context = {
        'results': results,
        'plot': plot,
        'problem': "-u'' + Pe·u' = 0,  u(0)=0, u(1)=1",
        'exact': 'u(x) = (exp(Pe·x) - 1) / (exp(Pe) - 1)',
    }
    return render(request, 'peclet_study.html', context)
