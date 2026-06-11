from django.shortcuts import render
from fenics import UnitIntervalMesh, UnitSquareMesh

from .forms import TransientForm
from .solver import (
    solve_transient, generate_transient_plotly, generate_transient_animation,
    generate_transient_snapshot_grid, scheme_comparison_1d, generate_scheme_comparison_plot,
)

SCHEME_LABELS = {
    'backward_euler': 'Backward Euler',
    'crank_nicolson': 'Crank-Nicolson',
}


def transient(request):
    if request.method == "POST":
        form = TransientForm(request.POST)
        if form.is_valid():
            d = form.cleaned_data
            dimension = d['dimension']
            n = d['n']

            mesh = UnitIntervalMesh(n) if dimension == '1D' else UnitSquareMesh(n, n)

            schemes = ['backward_euler', 'crank_nicolson'] if d['scheme'] == 'both' else [d['scheme']]
            results = []
            for scheme in schemes:
                times, snapshots, coords, _ = solve_transient(
                    dimension, mesh,
                    rho=d['rho'], cp=d['cp'], k=d['k'], f=d['f'],
                    T=d['T'], dt=d['dt'], scheme=scheme,
                    initial_condition=d['initial_condition'],
                    left_bc=d['left_bc'], right_bc=d['right_bc'],
                    bottom_bc=d['bottom_bc'], top_bc=d['top_bc'],
                )

                final_values = [[list(coords[i]), float(snapshots[-1][i])] for i in range(len(coords))]

                result = {
                    'scheme': scheme,
                    'label': SCHEME_LABELS[scheme],
                    'times': times,
                    'final_values': final_values,
                }

                if dimension == '1D':
                    result['plotly_json'] = generate_transient_plotly(times, snapshots, coords)
                    result['animation'] = generate_transient_animation(times, snapshots, coords, label=scheme)
                else:
                    result['snapshot_grid'] = generate_transient_snapshot_grid(
                        times, snapshots, coords, mesh.cells().tolist(), label=scheme,
                    )

                results.append(result)

            context = {
                'form': form,
                'results': results,
                'dimension': dimension,
                'nsteps': len(results[0]['times']) - 1,
                'message': 'Problem Solved',
            }
            return render(request, 'transient.html', context)
    else:
        form = TransientForm()

    return render(request, 'transient.html', {
        'form': form,
        'message': 'Enter Transient Heat Conduction Parameters',
    })


def scheme_comparison(request):
    """Temporal convergence study: Backward Euler vs Crank-Nicolson for the
    1D rod benchmark u(x,0)=sin(pi x), u(0,t)=u(1,t)=0, f=0, against the
    exact solution u(x,t) = sin(pi x) * exp(-pi^2 k t / (rho cp))."""
    results = scheme_comparison_1d()
    plot = generate_scheme_comparison_plot(results)

    be, cn = results['backward_euler'], results['crank_nicolson']
    rows = []
    for i, dt in enumerate(results['dt_values']):
        rows.append({
            'dt': dt,
            'be_error': be['errors'][i],
            'be_rate': be['rates'][i],
            'cn_error': cn['errors'][i],
            'cn_rate': cn['rates'][i],
        })

    context = {
        'rows': rows,
        'plot': plot,
        'problem': "rho*cp du/dt = k d²u/dx²,  u(x,0) = sin(pi x),  u(0,t) = u(1,t) = 0",
        'exact': 'u(x,t) = sin(pi x) · exp(-pi² k t / (rho cp))',
    }
    return render(request, 'scheme_comparison.html', context)
