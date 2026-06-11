"""Tests for the convection-diffusion solver — ConvectionDiffusion/solver.py

Covers Goal 2 (Convection-Diffusion Solver) and Goal 3 (SUPG Stabilization)
of the Phase 2 specification (docs/phases/phase2.md).
"""
import numpy as np
import pytest
from fenics import UnitIntervalMesh, UnitSquareMesh

from ConvectionDiffusion.solver import (
    compute_peclet_number,
    solve_convection_diffusion,
    mark_inflow_outflow,
    exact_1d_solution,
    peclet_study_1d,
    generate_convection_solution_plot,
    generate_peclet_plot,
)


# ---------------------------------------------------------------------------
# Cell Peclet number (Pe_h = |b| * h / (2k))
# ---------------------------------------------------------------------------

def test_peclet_number_diffusion_dominated():
    """k=1, b=1 on a 20-cell unit interval -> Pe_h < 1 (Galerkin stable)."""
    mesh = UnitIntervalMesh(20)
    pe = compute_peclet_number(mesh, k=1.0, bx=1.0)
    assert pe == pytest.approx(0.025, rel=1e-6)
    assert pe < 1


def test_peclet_number_convection_dominated():
    """k=0.001, b=1 on a 20-cell unit interval -> Pe_h > 1 (SUPG required)."""
    mesh = UnitIntervalMesh(20)
    pe = compute_peclet_number(mesh, k=0.001, bx=1.0)
    assert pe == pytest.approx(25.0, rel=1e-6)
    assert pe > 1


def test_peclet_number_zero_diffusivity_is_infinite():
    mesh = UnitIntervalMesh(20)
    assert compute_peclet_number(mesh, k=0.0, bx=1.0) == float('inf')


# ---------------------------------------------------------------------------
# Exact 1D solution: -u'' + Pe*u' = 0, u(0)=0, u(1)=1
# ---------------------------------------------------------------------------

def test_exact_1d_solution_satisfies_boundary_conditions():
    u_exact = exact_1d_solution(Pe=1.0)
    assert u_exact(0.0) == pytest.approx(0.0, abs=1e-12)
    assert u_exact(1.0) == pytest.approx(1.0, abs=1e-12)


def test_exact_1d_solution_zero_peclet_is_linear():
    """As Pe -> 0 the exact solution degenerates to u(x) = x."""
    u_exact = exact_1d_solution(Pe=0.0)
    assert u_exact(0.5) == pytest.approx(0.5, abs=1e-12)
    assert u_exact(0.25) == pytest.approx(0.25, abs=1e-12)


# ---------------------------------------------------------------------------
# Standard Galerkin correctness
# ---------------------------------------------------------------------------

def test_1d_galerkin_zero_velocity_reduces_to_pure_diffusion():
    """b=0 -> -k*u''=f, the same problem solved by SteadyStateThermal.

    f=2, k=1, u(0)=u(1)=0 is nodally exact for P1: u(x) = x(1-x).
    """
    mesh = UnitIntervalMesh(20)
    u_sol, _ = solve_convection_diffusion(
        '1D', mesh, k=1.0, bx=0.0, by=0.0, f=2.0,
        left_bc=0.0, right_bc=0.0, method='galerkin',
    )
    assert u_sol(0.5) == pytest.approx(0.25, abs=1e-10)


def test_1d_galerkin_low_peclet_matches_exact_solution():
    """At low cell Peclet number, standard Galerkin is close to the exact
    solution of -u'' + Pe*u' = 0."""
    mesh = UnitIntervalMesh(20)
    u_sol, _ = solve_convection_diffusion(
        '1D', mesh, k=1.0, bx=1.0, by=0.0, f=0.0,
        left_bc=0.0, right_bc=1.0, method='galerkin',
    )
    from fenics import errornorm
    u_exact = exact_1d_solution(Pe=1.0)
    err = errornorm(u_exact, u_sol, 'L2')
    assert err < 1e-3


# ---------------------------------------------------------------------------
# SUPG stabilization (Goal 3)
# ---------------------------------------------------------------------------

def test_1d_supg_reduces_to_pure_diffusion_when_velocity_is_zero():
    """The tau formula must not blow up (0/0) when b=0; the small-Pe limit
    tau -> h^2/(12k) recovers the same nodally exact pure-diffusion solution."""
    mesh = UnitIntervalMesh(20)
    u_sol, _ = solve_convection_diffusion(
        '1D', mesh, k=1.0, bx=0.0, by=0.0, f=2.0,
        left_bc=0.0, right_bc=0.0, method='supg',
    )
    assert u_sol(0.5) == pytest.approx(0.25, abs=1e-10)


def test_1d_galerkin_oscillates_at_high_cell_peclet():
    """At cell Pe_h=2.5 (k=0.01, b=1, h=0.05), unstabilized Galerkin must
    undershoot below the lower Dirichlet value (spurious oscillation)."""
    mesh = UnitIntervalMesh(20)
    u_sol, _ = solve_convection_diffusion(
        '1D', mesh, k=0.01, bx=1.0, by=0.0, f=0.0,
        left_bc=0.0, right_bc=1.0, method='galerkin',
    )
    vals = u_sol.compute_vertex_values(mesh)
    assert vals.min() < -0.1
    assert vals.max() == pytest.approx(1.0, abs=1e-8)


def test_1d_supg_remains_bounded_at_high_cell_peclet():
    """At the same cell Pe_h=2.5, SUPG-stabilized solution must respect the
    discrete maximum principle: stay within [0, 1]."""
    mesh = UnitIntervalMesh(20)
    u_sol, _ = solve_convection_diffusion(
        '1D', mesh, k=0.01, bx=1.0, by=0.0, f=0.0,
        left_bc=0.0, right_bc=1.0, method='supg',
    )
    vals = u_sol.compute_vertex_values(mesh)
    assert vals.min() >= -1e-8
    assert vals.max() == pytest.approx(1.0, abs=1e-8)


def test_1d_supg_more_accurate_than_galerkin_at_high_peclet():
    """At Pe=100 (cell Pe_h=2.5), SUPG's L2 error against the exact solution
    must be lower than standard Galerkin's (the convection-dominated regime
    where stabilization improves accuracy)."""
    results = peclet_study_1d([100])
    r = results[0]
    assert r['l2_supg'] < r['l2_galerkin']


# ---------------------------------------------------------------------------
# 2D solver
# ---------------------------------------------------------------------------

def test_2d_solver_returns_correct_structure():
    mesh = UnitSquareMesh(8, 8)
    u_sol, solution_list = solve_convection_diffusion(
        '2D', mesh, k=0.1, bx=1.0, by=0.0, f=0.0,
        left_bc=0.0, right_bc=1.0, bottom_bc=0.0, top_bc=0.0,
        method='galerkin',
    )
    assert u_sol is not None
    assert len(solution_list) == 81  # (8+1)^2 vertices
    first = solution_list[0]
    assert len(first) == 2
    assert len(first[0]) == 2
    assert isinstance(first[1], float)


def test_2d_solver_dirichlet_boundary_values():
    mesh = UnitSquareMesh(8, 8)
    u_sol, _ = solve_convection_diffusion(
        '2D', mesh, k=0.1, bx=1.0, by=0.0, f=0.0,
        left_bc=0.0, right_bc=1.0, bottom_bc=0.0, top_bc=0.0,
        method='galerkin',
    )
    assert u_sol(0.0, 0.5) == pytest.approx(0.0, abs=1e-8)
    assert u_sol(1.0, 0.5) == pytest.approx(1.0, abs=1e-8)


# ---------------------------------------------------------------------------
# Inflow/outflow boundary classification (Goal 2.6)
# ---------------------------------------------------------------------------

def test_mark_inflow_outflow_1d():
    """On the unit interval with b=1, the left endpoint (outward normal -1)
    is inflow (b.n<0) and the right endpoint (outward normal +1) is outflow."""
    mesh = UnitIntervalMesh(10)
    markers = mark_inflow_outflow(mesh, bx=1.0)
    counts = np.bincount(markers.array(), minlength=4)
    assert counts[1] == 1  # one inflow facet
    assert counts[2] == 1  # one outflow facet
    assert counts[0] == mesh.num_entities(0) - 2  # all other facets untouched


def test_mark_inflow_outflow_2d_with_velocity():
    """On a 4x4 unit square with b=(1,0): left edge (4 facets) is inflow,
    right edge (4 facets) is outflow, top+bottom (8 facets) are characteristic."""
    mesh = UnitSquareMesh(4, 4)
    markers = mark_inflow_outflow(mesh, bx=1.0, by=0.0)
    counts = np.bincount(markers.array(), minlength=4)
    assert counts[1] == 4
    assert counts[2] == 4
    assert counts[3] == 8


def test_mark_inflow_outflow_2d_zero_velocity_is_all_characteristic():
    """With b=0, b.n=0 everywhere -> every boundary facet is characteristic."""
    mesh = UnitSquareMesh(4, 4)
    markers = mark_inflow_outflow(mesh, bx=0.0, by=0.0)
    counts = np.bincount(markers.array(), minlength=4)
    assert counts[1] == 0
    assert counts[2] == 0
    assert counts[3] == 16


# ---------------------------------------------------------------------------
# Peclet number study (1D Galerkin vs. SUPG benchmark)
# ---------------------------------------------------------------------------

def test_peclet_study_1d_returns_one_result_per_pe():
    pe_values = [1, 10, 50, 100]
    results = peclet_study_1d(pe_values)
    assert len(results) == len(pe_values)
    for r, pe in zip(results, pe_values):
        assert r['pe'] == pe
        for key in ('cell_pe', 'x', 'u_exact', 'u_galerkin', 'u_supg',
                     'l2_galerkin', 'l2_supg'):
            assert key in r
        assert len(r['x']) == len(r['u_exact']) == len(r['u_galerkin']) == len(r['u_supg'])


def test_peclet_study_1d_supg_better_at_high_pe():
    results = peclet_study_1d([1, 10, 50, 100])
    high_pe = results[-1]
    assert high_pe['pe'] == 100
    assert high_pe['l2_supg'] < high_pe['l2_galerkin']


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def test_generate_convection_solution_plot_1d(tmp_path):
    mesh = UnitIntervalMesh(20)
    u_sol, _ = solve_convection_diffusion(
        '1D', mesh, k=0.01, bx=1.0, by=0.0, f=0.0,
        left_bc=0.0, right_bc=1.0, method='galerkin',
    )
    filename = generate_convection_solution_plot(
        u_sol, mesh, bx=1.0, by=0.0, label='galerkin', output_dir=str(tmp_path),
    )
    assert filename == 'convection_solution_galerkin.png'
    assert (tmp_path / filename).exists()


def test_generate_convection_solution_plot_2d_with_quiver(tmp_path):
    mesh = UnitSquareMesh(8, 8)
    u_sol, _ = solve_convection_diffusion(
        '2D', mesh, k=0.1, bx=1.0, by=1.0, f=0.0,
        left_bc=0.0, right_bc=1.0, bottom_bc=0.0, top_bc=0.0,
        method='supg',
    )
    filename = generate_convection_solution_plot(
        u_sol, mesh, bx=1.0, by=1.0, label='supg', output_dir=str(tmp_path),
    )
    assert filename == 'convection_solution_supg.png'
    assert (tmp_path / filename).exists()


def test_generate_convection_solution_plot_no_label(tmp_path):
    mesh = UnitIntervalMesh(10)
    u_sol, _ = solve_convection_diffusion(
        '1D', mesh, k=1.0, bx=0.0, by=0.0, f=0.0,
        left_bc=0.0, right_bc=1.0, method='galerkin',
    )
    filename = generate_convection_solution_plot(
        u_sol, mesh, bx=0.0, by=0.0, output_dir=str(tmp_path),
    )
    assert filename == 'convection_solution.png'
    assert (tmp_path / filename).exists()


def test_generate_peclet_plot(tmp_path):
    results = peclet_study_1d([1, 100])
    filename = generate_peclet_plot(results, output_dir=str(tmp_path))
    assert filename == 'peclet_study.png'
    assert (tmp_path / filename).exists()


# ---------------------------------------------------------------------------
# Django views
# ---------------------------------------------------------------------------

def test_convection_view_get_renders_form(client):
    response = client.get('/ConvectionDiffusion/')
    assert response.status_code == 200
    assert 'form' in response.context
    assert response.context['message'] == 'Enter Convection-Diffusion Parameters'


def test_convection_view_post_compares_galerkin_and_supg(client):
    data = {
        'dimension': '1D', 'n': 20,
        'k': 0.01, 'bx': 1.0, 'by': 0.0, 'f': 0.0,
        'left_bc': 0.0, 'right_bc': 1.0, 'bottom_bc': 0.0, 'top_bc': 0.0,
        'method': 'both',
    }
    response = client.post('/ConvectionDiffusion/', data)
    assert response.status_code == 200
    results = response.context['results']
    assert [r['method'] for r in results] == ['galerkin', 'supg']
    assert response.context['peclet'] == pytest.approx(2.5, rel=1e-6)


def test_convection_view_post_supg_only(client):
    data = {
        'dimension': '1D', 'n': 20,
        'k': 0.01, 'bx': 1.0, 'by': 0.0, 'f': 0.0,
        'left_bc': 0.0, 'right_bc': 1.0, 'bottom_bc': 0.0, 'top_bc': 0.0,
        'method': 'supg',
    }
    response = client.post('/ConvectionDiffusion/', data)
    assert response.status_code == 200
    results = response.context['results']
    assert len(results) == 1
    assert results[0]['method'] == 'supg'


def test_convection_view_post_2d(client):
    data = {
        'dimension': '2D', 'n': 8,
        'k': 0.1, 'bx': 1.0, 'by': 1.0, 'f': 0.0,
        'left_bc': 0.0, 'right_bc': 1.0, 'bottom_bc': 0.0, 'top_bc': 0.0,
        'method': 'galerkin',
    }
    response = client.post('/ConvectionDiffusion/', data)
    assert response.status_code == 200
    results = response.context['results']
    assert len(results) == 1
    assert len(results[0]['solution_list']) == 81


def test_peclet_study_view(client):
    response = client.get('/ConvectionDiffusion/peclet/')
    assert response.status_code == 200
    assert len(response.context['results']) == 4
    assert response.context['plot'] == 'peclet_study.png'
