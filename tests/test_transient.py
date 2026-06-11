"""Tests for the transient heat conduction solver — TransientThermal/solver.py

Covers Goal 4 (Transient Heat Conduction) of the Phase 3 specification
(docs/phases/phase3.md): rho*cp*du/dt = k*div(grad(u)) + f, integrated with
Backward Euler (1st-order) and Crank-Nicolson (2nd-order) time-stepping.
"""
import json

import numpy as np
import pytest
from fenics import UnitIntervalMesh, UnitSquareMesh, errornorm, Expression

from TransientThermal.solver import (
    solve_transient,
    exact_transient_1d,
    scheme_comparison_1d,
    generate_transient_plotly,
    generate_transient_animation,
    generate_transient_snapshot_grid,
    generate_scheme_comparison_plot,
)
from TransientThermal.forms import TransientForm, MAX_TIME_STEPS


# ---------------------------------------------------------------------------
# solve_transient — basic structure and the no-overshoot time-stepping fix
# ---------------------------------------------------------------------------

def test_solve_transient_1d_returns_correct_structure():
    mesh = UnitIntervalMesh(20)
    times, snapshots, coords, u_final = solve_transient(
        '1D', mesh, rho=1.0, cp=1.0, k=1.0, f=0.0, T=0.1, dt=0.02,
        scheme='backward_euler', initial_condition='sine',
        left_bc=0.0, right_bc=0.0,
    )
    assert times == pytest.approx([0.0, 0.02, 0.04, 0.06, 0.08, 0.1])
    assert len(snapshots) == len(times)
    assert len(coords) == 21  # n+1 vertices
    for snap in snapshots:
        assert len(snap) == 21
    assert u_final is not None


def test_solve_transient_nsteps_does_not_overshoot_T():
    """T/dt = 2.5 (non-integer) must produce nsteps = round(2.5) = 2 steps,
    i.e. times stop at 2*dt = 0.08, never overshooting T=0.1."""
    mesh = UnitIntervalMesh(20)
    times, snapshots, _, _ = solve_transient(
        '1D', mesh, rho=1.0, cp=1.0, k=1.0, f=0.0, T=0.1, dt=0.04,
        scheme='backward_euler', initial_condition='sine',
        left_bc=0.0, right_bc=0.0,
    )
    assert times == pytest.approx([0.0, 0.04, 0.08])
    assert all(t <= 0.1 for t in times)
    assert len(snapshots) == 3


def test_solve_transient_minimum_one_step():
    """Even if dt > T, at least one time step is taken."""
    mesh = UnitIntervalMesh(10)
    times, snapshots, _, _ = solve_transient(
        '1D', mesh, rho=1.0, cp=1.0, k=1.0, f=0.0, T=0.01, dt=0.02,
        scheme='backward_euler', initial_condition='sine',
        left_bc=0.0, right_bc=0.0,
    )
    assert len(times) == 2
    assert len(snapshots) == 2


# ---------------------------------------------------------------------------
# Steady-state preservation: if the initial condition already satisfies the
# steady-state equation (constant field matching the Dirichlet BCs, f=0),
# the solution must not change over time, for either scheme.
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('scheme', ['backward_euler', 'crank_nicolson'])
def test_constant_ic_remains_steady_state_1d(scheme):
    mesh = UnitIntervalMesh(20)
    _, snapshots, _, _ = solve_transient(
        '1D', mesh, rho=1.0, cp=1.0, k=1.0, f=0.0, T=0.1, dt=0.02,
        scheme=scheme, initial_condition='constant',
        left_bc=1.0, right_bc=1.0,
    )
    for snap in snapshots:
        assert np.allclose(snap, 1.0, atol=1e-9)


@pytest.mark.parametrize('scheme', ['backward_euler', 'crank_nicolson'])
def test_constant_ic_remains_steady_state_2d(scheme):
    mesh = UnitSquareMesh(6, 6)
    _, snapshots, _, _ = solve_transient(
        '2D', mesh, rho=1.0, cp=1.0, k=1.0, f=0.0, T=0.1, dt=0.02,
        scheme=scheme, initial_condition='constant',
        left_bc=1.0, right_bc=1.0, bottom_bc=1.0, top_bc=1.0,
    )
    for snap in snapshots:
        assert np.allclose(snap, 1.0, atol=1e-9)


# ---------------------------------------------------------------------------
# Initial conditions
# ---------------------------------------------------------------------------

def test_sine_ic_matches_exact_solution_at_t0():
    mesh = UnitIntervalMesh(40)
    times, snapshots, coords, _ = solve_transient(
        '1D', mesh, rho=1.0, cp=1.0, k=1.0, f=0.0, T=0.04, dt=0.02,
        scheme='backward_euler', initial_condition='sine',
        left_bc=0.0, right_bc=0.0,
    )
    coords = np.array(coords).flatten()
    expected = np.sin(np.pi * coords)
    assert np.allclose(snapshots[0], expected, atol=1e-12)


def test_gaussian_ic_is_bounded_in_unit_interval():
    """The Gaussian pulse exp(-50*(x-0.5)^2) is bounded in (0, 1] for x in
    [0, 1], and diffusion can only spread/decay it. P1 Galerkin does not
    strictly satisfy the discrete maximum principle on this 8x8 mesh for
    such a steep profile, producing a small Gibbs-type undershoot, so allow
    a small tolerance below 0 rather than requiring strict boundedness."""
    mesh = UnitSquareMesh(8, 8)
    times, snapshots, _, _ = solve_transient(
        '2D', mesh, rho=1.0, cp=1.0, k=1.0, f=0.0, T=0.05, dt=0.01,
        scheme='crank_nicolson', initial_condition='gaussian',
        left_bc=0.0, right_bc=0.0, bottom_bc=0.0, top_bc=0.0,
    )
    snapshots = np.array(snapshots)
    assert snapshots.min() >= -0.01
    assert snapshots.max() <= 1.0 + 1e-8
    assert snapshots[0].max() == pytest.approx(1.0, abs=1e-6)


# ---------------------------------------------------------------------------
# Accuracy: Backward Euler (O(dt)) vs Crank-Nicolson (O(dt^2)) against the
# analytical decay solution u(x,t) = sin(pi x) * exp(-pi^2 k t / (rho cp))
# ---------------------------------------------------------------------------

def test_backward_euler_approximates_exact_decay():
    mesh = UnitIntervalMesh(50)
    _, _, _, u_final = solve_transient(
        '1D', mesh, rho=1.0, cp=1.0, k=1.0, f=0.0, T=0.1, dt=0.01,
        scheme='backward_euler', initial_condition='sine',
        left_bc=0.0, right_bc=0.0,
    )
    u_exact_T = Expression(
        "sin(pi*x[0]) * exp(-pi*pi*k*t/(rho*cp))", degree=4,
        k=1.0, rho=1.0, cp=1.0, t=0.1,
    )
    err = errornorm(u_exact_T, u_final, 'L2')
    assert err < 0.02


def test_crank_nicolson_more_accurate_than_backward_euler_at_same_dt():
    mesh = UnitIntervalMesh(50)
    u_exact_T = Expression(
        "sin(pi*x[0]) * exp(-pi*pi*k*t/(rho*cp))", degree=4,
        k=1.0, rho=1.0, cp=1.0, t=0.1,
    )

    _, _, _, u_be = solve_transient(
        '1D', mesh, rho=1.0, cp=1.0, k=1.0, f=0.0, T=0.1, dt=0.01,
        scheme='backward_euler', initial_condition='sine',
        left_bc=0.0, right_bc=0.0,
    )
    _, _, _, u_cn = solve_transient(
        '1D', mesh, rho=1.0, cp=1.0, k=1.0, f=0.0, T=0.1, dt=0.01,
        scheme='crank_nicolson', initial_condition='sine',
        left_bc=0.0, right_bc=0.0,
    )

    err_be = errornorm(u_exact_T, u_be, 'L2')
    err_cn = errornorm(u_exact_T, u_cn, 'L2')
    assert err_cn < err_be / 5


# ---------------------------------------------------------------------------
# exact_transient_1d
# ---------------------------------------------------------------------------

def test_exact_transient_1d_at_t0_is_initial_condition():
    x = np.linspace(0, 1, 11)
    assert np.allclose(exact_transient_1d(x, 0.0), np.sin(np.pi * x))


def test_exact_transient_1d_decays_over_time():
    u0 = exact_transient_1d(0.5, 0.0)
    u1 = exact_transient_1d(0.5, 0.1)
    assert 0 < u1 < u0


# ---------------------------------------------------------------------------
# scheme_comparison_1d — temporal convergence study
# ---------------------------------------------------------------------------

def test_scheme_comparison_1d_structure():
    results = scheme_comparison_1d(dt_values=(0.02, 0.01, 0.005), T=0.1, n=50)
    assert results['dt_values'] == [0.02, 0.01, 0.005]
    for scheme in ('backward_euler', 'crank_nicolson'):
        assert len(results[scheme]['errors']) == 3
        assert len(results[scheme]['rates']) == 3
        assert results[scheme]['rates'][0] is None
        assert len(results[scheme]['times']) == len(results[scheme]['u_center'])
    assert len(results['exact_u_center']) == len(results['backward_euler']['times'])


def test_scheme_comparison_1d_backward_euler_rate_near_one():
    results = scheme_comparison_1d(dt_values=(0.02, 0.01, 0.005), T=0.1, n=200)
    rates = results['backward_euler']['rates'][1:]
    for rate in rates:
        assert 0.8 < rate < 1.2


def test_scheme_comparison_1d_crank_nicolson_rate_near_two():
    results = scheme_comparison_1d(dt_values=(0.02, 0.01, 0.005), T=0.1, n=200)
    rates = results['crank_nicolson']['rates'][1:]
    for rate in rates:
        assert 1.5 < rate < 2.1


def test_scheme_comparison_1d_crank_nicolson_more_accurate():
    results = scheme_comparison_1d(dt_values=(0.02, 0.01, 0.005), T=0.1, n=200)
    for be_err, cn_err in zip(results['backward_euler']['errors'], results['crank_nicolson']['errors']):
        assert cn_err < be_err


# ---------------------------------------------------------------------------
# Plotting / export helpers
# ---------------------------------------------------------------------------

def test_generate_transient_plotly_returns_valid_json_with_frames():
    mesh = UnitIntervalMesh(20)
    times, snapshots, coords, _ = solve_transient(
        '1D', mesh, rho=1.0, cp=1.0, k=1.0, f=0.0, T=0.04, dt=0.02,
        scheme='backward_euler', initial_condition='sine',
        left_bc=0.0, right_bc=0.0,
    )
    fig_json = generate_transient_plotly(times, snapshots, coords)
    fig = json.loads(fig_json)
    assert 'data' in fig and 'layout' in fig and 'frames' in fig
    assert len(fig['frames']) == len(times)
    assert 'sliders' in fig['layout']


def test_generate_transient_animation_creates_gif(tmp_path):
    mesh = UnitIntervalMesh(20)
    times, snapshots, coords, _ = solve_transient(
        '1D', mesh, rho=1.0, cp=1.0, k=1.0, f=0.0, T=0.04, dt=0.02,
        scheme='backward_euler', initial_condition='sine',
        left_bc=0.0, right_bc=0.0,
    )
    filename = generate_transient_animation(
        times, snapshots, coords, label='backward_euler', output_dir=str(tmp_path),
    )
    assert filename == 'transient_animation_backward_euler.gif'
    assert (tmp_path / filename).exists()


def test_generate_transient_snapshot_grid_creates_png(tmp_path):
    mesh = UnitSquareMesh(6, 6)
    times, snapshots, coords, _ = solve_transient(
        '2D', mesh, rho=1.0, cp=1.0, k=1.0, f=0.0, T=0.04, dt=0.02,
        scheme='crank_nicolson', initial_condition='gaussian',
        left_bc=0.0, right_bc=0.0, bottom_bc=0.0, top_bc=0.0,
    )
    filename = generate_transient_snapshot_grid(
        times, snapshots, coords, mesh.cells().tolist(),
        label='crank_nicolson', output_dir=str(tmp_path),
    )
    assert filename == 'transient_snapshots_crank_nicolson.png'
    assert (tmp_path / filename).exists()


def test_generate_scheme_comparison_plot_creates_png(tmp_path):
    results = scheme_comparison_1d(dt_values=(0.02, 0.01), T=0.1, n=50)
    filename = generate_scheme_comparison_plot(results, output_dir=str(tmp_path))
    assert filename == 'scheme_comparison.png'
    assert (tmp_path / filename).exists()


# ---------------------------------------------------------------------------
# Forms — resource guards (Goal 4: timestep controls)
# ---------------------------------------------------------------------------

def _form_data(**overrides):
    data = {
        'dimension': '1D', 'n': 20,
        'rho': 1.0, 'cp': 1.0, 'k': 1.0, 'f': 0.0,
        'T': 0.1, 'dt': 0.01,
        'initial_condition': 'sine', 'scheme': 'both',
        'left_bc': 0.0, 'right_bc': 0.0, 'bottom_bc': 0.0, 'top_bc': 0.0,
    }
    data.update(overrides)
    return data


def test_form_valid_with_defaults():
    form = TransientForm(_form_data())
    assert form.is_valid(), form.errors


def test_form_rejects_dt_greater_than_T():
    form = TransientForm(_form_data(T=0.01, dt=0.1))
    assert not form.is_valid()


def test_form_rejects_too_many_time_steps():
    form = TransientForm(_form_data(T=10.0, dt=0.001))
    assert round(10.0 / 0.001) > MAX_TIME_STEPS
    assert not form.is_valid()


def test_form_defaults_2d_only_fields_when_blank():
    data = _form_data(dimension='2D')
    del data['bottom_bc']
    del data['top_bc']
    form = TransientForm(data)
    assert form.is_valid(), form.errors
    assert form.cleaned_data['bottom_bc'] == 0.0
    assert form.cleaned_data['top_bc'] == 0.0


# ---------------------------------------------------------------------------
# Django views
# ---------------------------------------------------------------------------

def test_transient_view_get_renders_form(client):
    response = client.get('/TransientThermal/')
    assert response.status_code == 200
    assert 'form' in response.context
    assert response.context['message'] == 'Enter Transient Heat Conduction Parameters'


def test_transient_view_post_1d_both_schemes(client):
    response = client.post('/TransientThermal/', _form_data())
    assert response.status_code == 200
    results = response.context['results']
    assert [r['scheme'] for r in results] == ['backward_euler', 'crank_nicolson']
    for r in results:
        assert 'plotly_json' in r
        assert 'animation' in r
        assert len(r['final_values']) == 21
    assert response.context['nsteps'] == 10


def test_transient_view_post_1d_single_scheme(client):
    response = client.post('/TransientThermal/', _form_data(scheme='crank_nicolson'))
    assert response.status_code == 200
    results = response.context['results']
    assert len(results) == 1
    assert results[0]['scheme'] == 'crank_nicolson'


def test_transient_view_post_2d(client):
    data = _form_data(dimension='2D', n=8, scheme='backward_euler',
                       left_bc=0.0, right_bc=1.0, bottom_bc=0.0, top_bc=0.0)
    response = client.post('/TransientThermal/', data)
    assert response.status_code == 200
    results = response.context['results']
    assert len(results) == 1
    assert 'snapshot_grid' in results[0]
    assert len(results[0]['final_values']) == 81  # (8+1)^2 vertices


def test_transient_view_post_invalid_form_reshows_form(client):
    response = client.post('/TransientThermal/', _form_data(T=0.01, dt=0.1))
    assert response.status_code == 200
    assert not response.context['form'].is_valid()
    assert 'results' not in response.context


def test_scheme_comparison_view(client):
    response = client.get('/TransientThermal/comparison/')
    assert response.status_code == 200
    assert response.context['plot'] == 'scheme_comparison.png'
    rows = response.context['rows']
    assert len(rows) == 3
    assert rows[0]['be_rate'] is None
    assert rows[1]['be_rate'] is not None
