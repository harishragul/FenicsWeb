"""Mesh-refinement (h) and polynomial-degree-refinement (p) convergence
studies for the Validation module.

Implements Phase 4 tasks 4.2-4.4 of docs/phases/phase4.md:

- 1D and 2D h-refinement for the Poisson equation (P1), reporting L2 and H1
  errors and observed convergence rates against the manufactured solutions
  in `Validation.analytical`.
- p-refinement of the 2D manufactured-solution problem with P1/P2/P3
  elements.
- Galerkin vs. SUPG convergence for the 1D convection-diffusion benchmark at
  low (Pe=0.1) and high (Pe=100) Peclet number, reusing
  `ConvectionDiffusion.solver`.

All convergence rates are computed as

    rate = log(e_coarse / e_fine) / log(h_coarse / h_fine)

per docs/science/governing_equations.md Section 4.
"""
import math
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from fenics import (
    Constant, DirichletBC, Function, FunctionSpace, TestFunction,
    TrialFunction, UnitIntervalMesh, UnitSquareMesh, dot, dx, errornorm,
    grad, solve,
)

from FenicsWeb.settings import BASE_DIR
from ConvectionDiffusion.solver import (
    compute_peclet_number, exact_1d_solution, solve_convection_diffusion,
)
from .analytical import exact_1d_poisson, exact_2d_poisson, source_2d_poisson

STATIC_DIR = os.path.join(BASE_DIR, 'static')


def _rate(prev_err, err, prev_h, h):
    """Observed convergence rate between two successive mesh refinements."""
    if prev_err is None or err <= 0 or prev_err <= 0:
        return None
    return math.log(prev_err / err) / math.log(prev_h / h)


def poisson_1d_convergence(mesh_sizes=(4, 8, 16, 32, 64), k=1.0):
    """h-refinement study for -k*u'' = 2k on [0,1], u(0)=u(1)=0.

    Exact solution: u(x) = x(1-x). Returns a list of dicts with L2/H1 errors
    and observed rates (P1 elements).
    """
    u_exact = exact_1d_poisson()
    rows = []
    prev_l2 = prev_h1 = prev_h = None
    for n in mesh_sizes:
        mesh = UnitIntervalMesh(n)
        V = FunctionSpace(mesh, "P", 1)
        u = TrialFunction(V)
        v = TestFunction(V)
        bcs = [DirichletBC(V, Constant(0.0), "on_boundary")]
        a = k * dot(grad(u), grad(v)) * dx
        L = Constant(2.0 * k) * v * dx
        u_sol = Function(V)
        solve(a == L, u_sol, bcs)

        h = 1.0 / n
        l2 = errornorm(u_exact, u_sol, 'L2')
        h1 = errornorm(u_exact, u_sol, 'H1')
        rows.append({
            'n': n,
            'h': h,
            'l2_error': l2,
            'l2_rate': _rate(prev_l2, l2, prev_h, h),
            'h1_error': h1,
            'h1_rate': _rate(prev_h1, h1, prev_h, h),
        })
        prev_l2, prev_h1, prev_h = l2, h1, h
    return rows


def poisson_2d_convergence(mesh_sizes=(4, 8, 16, 32, 64), k=1.0, element_order=1):
    """h-refinement study for -k*lap(u) = f on the unit square,
    manufactured solution u(x,y) = sin(pi x) sin(pi y), u=0 on the boundary.

    Returns a list of dicts with L2/H1 errors and observed rates for the
    given element order (P1, P2, or P3).
    """
    u_exact = exact_2d_poisson()
    f_expr = source_2d_poisson(k=k)
    rows = []
    prev_l2 = prev_h1 = prev_h = None
    for n in mesh_sizes:
        mesh = UnitSquareMesh(n, n)
        V = FunctionSpace(mesh, "P", element_order)
        u = TrialFunction(V)
        v = TestFunction(V)
        bcs = [DirichletBC(V, Constant(0.0), "on_boundary")]
        a = k * dot(grad(u), grad(v)) * dx
        L = f_expr * v * dx
        u_sol = Function(V)
        solve(a == L, u_sol, bcs)

        h = 1.0 / n
        l2 = errornorm(u_exact, u_sol, 'L2')
        h1 = errornorm(u_exact, u_sol, 'H1')
        rows.append({
            'n': n,
            'h': h,
            'l2_error': l2,
            'l2_rate': _rate(prev_l2, l2, prev_h, h),
            'h1_error': h1,
            'h1_rate': _rate(prev_h1, h1, prev_h, h),
        })
        prev_l2, prev_h1, prev_h = l2, h1, h
    return rows


def p_refinement_study(mesh_sizes=(4, 8, 16, 32), orders=(1, 2, 3), k=1.0):
    """p-refinement study: solve the 2D manufactured-solution problem with
    P1, P2, and P3 elements across `mesh_sizes`.

    Returns a dict {element_order: rows}, where `rows` has the same
    structure as `poisson_2d_convergence`. Expected L2 rates are
    P1≈2, P2≈3, P3≈4 (docs/science/benchmark_references.md Section 3).
    """
    return {p: poisson_2d_convergence(mesh_sizes, k=k, element_order=p) for p in orders}


def convection_diffusion_convergence(pe_values=(0.1, 100), mesh_sizes=(8, 16, 32, 64)):
    """Galerkin vs. SUPG convergence for the 1D convection-diffusion
    benchmark -u'' + Pe*u' = 0, u(0)=0, u(1)=1.

    For each Peclet number and mesh size, solves with both `'galerkin'` and
    `'supg'` and reports the L2 error against the exact solution as well as
    the min/max nodal values (used to detect spurious over-/under-shoots
    when the cell Peclet number Pe_h > 1).

    Returns a dict {pe: rows}.
    """
    results = {}
    for pe in pe_values:
        k = 1.0 / pe
        u_exact = exact_1d_solution(pe)
        rows = []
        for n in mesh_sizes:
            mesh = UnitIntervalMesh(n)
            pe_h = compute_peclet_number(mesh, k, 1.0, 0.0)
            row = {'n': n, 'h': 1.0 / n, 'pe_h': pe_h}
            for method in ('galerkin', 'supg'):
                u_sol, _ = solve_convection_diffusion(
                    '1D', mesh, k=k, bx=1.0, by=0.0, f=0.0,
                    left_bc=0.0, right_bc=1.0, method=method,
                )
                vals = u_sol.compute_vertex_values(mesh)
                row[f'l2_{method}'] = errornorm(u_exact, u_sol, 'L2')
                row[f'min_{method}'] = float(vals.min())
                row[f'max_{method}'] = float(vals.max())
            rows.append(row)
        results[pe] = rows
    return results


def _plot_series(ax, h_values, err_values, label, color, marker, ref_slope):
    """Plot one error series on log-log axes with a dashed reference slope
    line and an annotation of the observed convergence rate."""
    ax.loglog(h_values, err_values, marker=marker, color=color, label=label)

    observed = _rate(err_values[0], err_values[-1], h_values[0], h_values[-1])
    if observed is not None:
        ax.lines[-1].set_label(f"{label} (observed ≈ {observed:.2f})")

    if ref_slope is not None:
        h_ref = np.array([h_values[0], h_values[-1]])
        c = err_values[0] / (h_values[0] ** ref_slope)
        ref = c * h_ref ** ref_slope
        ax.loglog(h_ref, ref, '--', color=color, alpha=0.35, linewidth=1)
        ax.annotate(f"O(h^{ref_slope})", xy=(h_ref[-1], ref[-1]),
                     fontsize=8, color=color, alpha=0.8,
                     xytext=(3, 0), textcoords='offset points')


def generate_h_refinement_plot(rows_1d, rows_2d, output_dir=None):
    """Log-log convergence plot (error vs. h) for the 1D and 2D Poisson
    h-refinement studies, with O(h) reference slopes annotated."""
    target_dir = output_dir if output_dir is not None else STATIC_DIR
    os.makedirs(target_dir, exist_ok=True)
    plot_filename = os.path.join(target_dir, 'h_refinement.png')

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    h1d = [r['h'] for r in rows_1d]
    _plot_series(axes[0], h1d, [r['l2_error'] for r in rows_1d], 'L2', 'C0', 'o', 2)
    _plot_series(axes[0], h1d, [r['h1_error'] for r in rows_1d], 'H1', 'C1', 's', 1)
    axes[0].set_title('1D Poisson h-refinement (P1)')
    axes[0].set_xlabel('h')
    axes[0].set_ylabel('error')
    axes[0].legend(fontsize=8)
    axes[0].grid(True, which='both', alpha=0.3)

    h2d = [r['h'] for r in rows_2d]
    _plot_series(axes[1], h2d, [r['l2_error'] for r in rows_2d], 'L2', 'C0', 'o', 2)
    _plot_series(axes[1], h2d, [r['h1_error'] for r in rows_2d], 'H1', 'C1', 's', 1)
    axes[1].set_title('2D Poisson h-refinement (P1)')
    axes[1].set_xlabel('h')
    axes[1].set_ylabel('error')
    axes[1].legend(fontsize=8)
    axes[1].grid(True, which='both', alpha=0.3)

    plt.tight_layout()
    plt.savefig(plot_filename, dpi=120)
    plt.close(fig)
    return 'h_refinement.png'


def generate_p_refinement_plot(p_results, output_dir=None):
    """Log-log L2-error convergence plot for the p-refinement study
    (P1/P2/P3), with O(h^{p+1}) reference slopes annotated."""
    target_dir = output_dir if output_dir is not None else STATIC_DIR
    os.makedirs(target_dir, exist_ok=True)
    plot_filename = os.path.join(target_dir, 'p_refinement.png')

    fig, ax = plt.subplots(figsize=(6, 4.5))
    colors = {1: 'C0', 2: 'C1', 3: 'C2'}
    markers = {1: 'o', 2: 's', 3: '^'}
    for order, rows in sorted(p_results.items()):
        h = [r['h'] for r in rows]
        l2 = [r['l2_error'] for r in rows]
        _plot_series(ax, h, l2, f'P{order}', colors.get(order, 'C3'),
                      markers.get(order, 'd'), order + 1)

    ax.set_title('p-refinement: 2D manufactured solution')
    ax.set_xlabel('h')
    ax.set_ylabel('L2 error')
    ax.legend(fontsize=8)
    ax.grid(True, which='both', alpha=0.3)

    plt.tight_layout()
    plt.savefig(plot_filename, dpi=120)
    plt.close(fig)
    return 'p_refinement.png'


def generate_convection_diffusion_plot(cd_results, output_dir=None):
    """Log-log L2-error convergence plot (Galerkin vs. SUPG) for each
    Peclet number in `cd_results`."""
    target_dir = output_dir if output_dir is not None else STATIC_DIR
    os.makedirs(target_dir, exist_ok=True)
    plot_filename = os.path.join(target_dir, 'cd_convergence.png')

    pe_values = sorted(cd_results.keys())
    fig, axes = plt.subplots(1, len(pe_values), figsize=(5.5 * len(pe_values), 4.5))
    if len(pe_values) == 1:
        axes = [axes]

    for ax, pe in zip(axes, pe_values):
        rows = cd_results[pe]
        h = [r['h'] for r in rows]
        _plot_series(ax, h, [r['l2_galerkin'] for r in rows], 'Galerkin', 'C3', 'o', None)
        _plot_series(ax, h, [r['l2_supg'] for r in rows], 'SUPG', 'C0', 's', None)
        ax.set_title(f'1D convection-diffusion, Pe = {pe}')
        ax.set_xlabel('h')
        ax.set_ylabel('L2 error')
        ax.legend(fontsize=8)
        ax.grid(True, which='both', alpha=0.3)

    plt.tight_layout()
    plt.savefig(plot_filename, dpi=120)
    plt.close(fig)
    return 'cd_convergence.png'
