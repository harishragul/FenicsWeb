from django.shortcuts import render

from .benchmarks import (
    REFERENCE_NUSSELT, de_vahl_davis_conduction, generate_de_vahl_davis_plot,
)
from .convergence import (
    convection_diffusion_convergence, generate_convection_diffusion_plot,
    generate_h_refinement_plot, generate_p_refinement_plot,
    p_refinement_study, poisson_1d_convergence, poisson_2d_convergence,
)


def _check(name, expected, observed, passed):
    return {'name': name, 'expected': expected, 'observed': observed, 'passed': passed}


def _fmt_rate(x):
    return f'{x:.2f}' if x is not None else '—'


def _fmt_rows(rows):
    """Format h/p-refinement rows for display (errors in scientific
    notation, which Django's `floatformat` filter cannot produce)."""
    return [
        {
            'n': r['n'],
            'h': f"{r['h']:.4f}",
            'l2_error': f"{r['l2_error']:.3e}",
            'l2_rate': _fmt_rate(r['l2_rate']),
            'h1_error': f"{r['h1_error']:.3e}",
            'h1_rate': _fmt_rate(r['h1_rate']),
        }
        for r in rows
    ]


def _fmt_cd_rows(rows):
    """Format convection-diffusion convergence rows for display."""
    return [
        {
            'n': r['n'],
            'h': f"{r['h']:.4f}",
            'pe_h': f"{r['pe_h']:.3f}",
            'l2_galerkin': f"{r['l2_galerkin']:.3e}",
            'l2_supg': f"{r['l2_supg']:.3e}",
            'min_galerkin': f"{r['min_galerkin']:.3f}",
            'max_galerkin': f"{r['max_galerkin']:.3f}",
            'min_supg': f"{r['min_supg']:.3f}",
            'max_supg': f"{r['max_supg']:.3f}",
        }
        for r in rows
    ]


def convergence(request):
    """Convergence (verification) studies: h-refinement (1D & 2D), p-refinement,
    and Galerkin vs. SUPG convergence for the convection-diffusion benchmark.
    """
    rows_1d = poisson_1d_convergence()
    rows_2d = poisson_2d_convergence()
    p_results = p_refinement_study()
    cd_results = convection_diffusion_convergence()

    context = {
        'rows_1d': _fmt_rows(rows_1d),
        'rows_2d': _fmt_rows(rows_2d),
        'p_results': [(order, _fmt_rows(rows)) for order, rows in sorted(p_results.items())],
        'cd_results': [(pe, _fmt_cd_rows(rows)) for pe, rows in sorted(cd_results.items())],
        'plot_h': generate_h_refinement_plot(rows_1d, rows_2d),
        'plot_p': generate_p_refinement_plot(p_results),
        'plot_cd': generate_convection_diffusion_plot(cd_results),
    }
    return render(request, 'convergence.html', context)


def benchmark(request):
    """De Vahl Davis cavity benchmark (pure-conduction / Ra=0 limit) with
    average Nusselt number compared to the reference table."""
    u_sol, mesh, nu_avg = de_vahl_davis_conduction()
    context = {
        'nu_avg': nu_avg,
        'reference': REFERENCE_NUSSELT,
        'plot': generate_de_vahl_davis_plot(u_sol, mesh),
    }
    return render(request, 'benchmark.html', context)


def dashboard(request):
    """Validation dashboard: a single page summarizing every verification
    and benchmark check with pass/fail indicators against the expected
    theoretical results documented in docs/science/benchmark_references.md.
    """
    rows_1d = poisson_1d_convergence(mesh_sizes=(8, 16, 32))
    rows_2d = poisson_2d_convergence(mesh_sizes=(8, 16, 32))
    p_results = p_refinement_study(mesh_sizes=(8, 16, 32))
    cd_results = convection_diffusion_convergence(mesh_sizes=(8, 16, 32))
    _, _, nu_avg = de_vahl_davis_conduction(n=16)

    l2_1d, h1_1d = rows_1d[-1]['l2_rate'], rows_1d[-1]['h1_rate']
    l2_2d, h1_2d = rows_2d[-1]['l2_rate'], rows_2d[-1]['h1_rate']
    p2_rate = p_results[2][-1]['l2_rate']
    p3_rate = p_results[3][-1]['l2_rate']
    cd100 = cd_results[100][0]  # coarsest mesh: Pe_h > 1, Galerkin oscillates

    checks = [
        _check('1D Poisson — L2 rate (P1)', '≈ 2.0', f'{l2_1d:.2f}',
               abs(l2_1d - 2.0) < 0.15),
        _check('1D Poisson — H1 rate (P1)', '≈ 1.0', f'{h1_1d:.2f}',
               abs(h1_1d - 1.0) < 0.1),
        _check('2D Poisson — L2 rate (P1)', '≈ 2.0', f'{l2_2d:.2f}',
               abs(l2_2d - 2.0) < 0.15),
        _check('2D Poisson — H1 rate (P1)', '≈ 1.0', f'{h1_2d:.2f}',
               abs(h1_2d - 1.0) < 0.1),
        _check('p-refinement — P2 L2 rate', '≈ 3.0', f'{p2_rate:.2f}',
               abs(p2_rate - 3.0) < 0.25),
        _check('p-refinement — P3 L2 rate', '≈ 4.0', f'{p3_rate:.2f}',
               abs(p3_rate - 4.0) < 0.4),
        _check(
            'Convection-diffusion — SUPG bounded at Pe=100, Pe_h > 1',
            'SUPG ∈ [0,1], Galerkin under/overshoots',
            f"Galerkin min={cd100['min_galerkin']:.3f}, SUPG min={cd100['min_supg']:.3f}",
            cd100['min_supg'] > -0.05 and cd100['min_galerkin'] < -0.05,
        ),
        _check('De Vahl Davis — Nu (Ra=0, pure conduction)', '= 1.0',
               f'{nu_avg:.4f}', abs(nu_avg - 1.0) < 0.01),
    ]
    return render(request, 'dashboard.html', {'checks': checks})
