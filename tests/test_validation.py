"""Tests for the Validation module (Phase 4 — Verification Framework and
Benchmark Library, docs/phases/phase4.md, Goals 6 & 7 of the project
specification).

Covers:
- analytical exact-solution expressions
- h-refinement convergence (1D & 2D Poisson, P1): L2 -> O(h^2), H1 -> O(h^1)
- p-refinement convergence (P1/P2/P3): L2 -> O(h^{p+1})
- Galerkin vs. SUPG convergence for the 1D convection-diffusion benchmark
- the De Vahl Davis pure-conduction (Ra=0) Nusselt-number benchmark
- plotting helpers
- Django views (dashboard, convergence, benchmark)
"""
import math

import pytest
from fenics import UnitIntervalMesh, UnitSquareMesh

from Validation.analytical import (
    exact_1d_poisson, exact_2d_poisson, source_2d_poisson,
)
from Validation.convergence import (
    _rate,
    convection_diffusion_convergence,
    generate_convection_diffusion_plot,
    generate_h_refinement_plot,
    generate_p_refinement_plot,
    p_refinement_study,
    poisson_1d_convergence,
    poisson_2d_convergence,
)
from Validation.benchmarks import (
    REFERENCE_NUSSELT,
    de_vahl_davis_conduction,
    generate_de_vahl_davis_plot,
)
from Validation.views import _check, _fmt_cd_rows, _fmt_rate, _fmt_rows


# ---------------------------------------------------------------------------
# analytical.py — exact (manufactured) solutions
# ---------------------------------------------------------------------------

def test_exact_1d_poisson_values():
    u = exact_1d_poisson()
    assert u(0.0) == pytest.approx(0.0)
    assert u(1.0) == pytest.approx(0.0)
    assert u(0.5) == pytest.approx(0.25)
    assert u(0.25) == pytest.approx(0.25 * 0.75)


def test_exact_2d_poisson_values():
    u = exact_2d_poisson()
    assert u(0.5, 0.5) == pytest.approx(1.0)
    assert u(0.0, 0.5) == pytest.approx(0.0)
    assert u(0.5, 0.0) == pytest.approx(0.0)
    assert u(1.0, 1.0) == pytest.approx(0.0, abs=1e-12)


def test_source_2d_poisson_values():
    f = source_2d_poisson(k=1.0)
    assert f(0.5, 0.5) == pytest.approx(2 * math.pi**2)

    f2 = source_2d_poisson(k=2.0)
    assert f2(0.5, 0.5) == pytest.approx(2 * 2 * math.pi**2)


# ---------------------------------------------------------------------------
# convergence.py — _rate helper
# ---------------------------------------------------------------------------

def test_rate_returns_none_for_first_point():
    assert _rate(None, 1e-3, 0.25, 0.125) is None


def test_rate_computes_log_ratio():
    # error halves when h halves -> rate = 1
    assert _rate(2e-3, 1e-3, 0.25, 0.125) == pytest.approx(1.0)
    # error quarters when h halves -> rate = 2
    assert _rate(4e-3, 1e-3, 0.25, 0.125) == pytest.approx(2.0)


def test_rate_returns_none_for_nonpositive_error():
    assert _rate(0.0, 1e-3, 0.25, 0.125) is None
    assert _rate(1e-3, 0.0, 0.25, 0.125) is None


# ---------------------------------------------------------------------------
# 1D Poisson h-refinement (Task 4.2)
# ---------------------------------------------------------------------------

def test_poisson_1d_convergence_structure():
    rows = poisson_1d_convergence(mesh_sizes=(4, 8, 16))
    assert len(rows) == 3
    assert rows[0]['l2_rate'] is None
    assert rows[0]['h1_rate'] is None
    for r in rows:
        assert set(r.keys()) == {'n', 'h', 'l2_error', 'l2_rate', 'h1_error', 'h1_rate'}
        assert r['h'] == pytest.approx(1.0 / r['n'])


def test_poisson_1d_convergence_rates_match_p1_theory():
    rows = poisson_1d_convergence(mesh_sizes=(8, 16, 32, 64))
    # L2 -> O(h^2), H1 -> O(h^1) for P1 elements
    assert rows[-1]['l2_rate'] == pytest.approx(2.0, abs=0.05)
    assert rows[-1]['h1_rate'] == pytest.approx(1.0, abs=0.05)
    # errors should also be small (close to the exact solution)
    assert rows[-1]['l2_error'] < 1e-3


# ---------------------------------------------------------------------------
# 2D Poisson h-refinement (Task 4.2)
# ---------------------------------------------------------------------------

def test_poisson_2d_convergence_structure():
    rows = poisson_2d_convergence(mesh_sizes=(4, 8, 16))
    assert len(rows) == 3
    assert rows[0]['l2_rate'] is None
    for r in rows:
        assert set(r.keys()) == {'n', 'h', 'l2_error', 'l2_rate', 'h1_error', 'h1_rate'}


def test_poisson_2d_convergence_rates_match_p1_theory():
    rows = poisson_2d_convergence(mesh_sizes=(8, 16, 32))
    assert rows[-1]['l2_rate'] == pytest.approx(2.0, abs=0.1)
    assert rows[-1]['h1_rate'] == pytest.approx(1.0, abs=0.05)


# ---------------------------------------------------------------------------
# p-refinement study (Task 4.3, Goal 5)
# ---------------------------------------------------------------------------

def test_p_refinement_study_structure():
    results = p_refinement_study(mesh_sizes=(4, 8, 16), orders=(1, 2, 3))
    assert set(results.keys()) == {1, 2, 3}
    for order, rows in results.items():
        assert len(rows) == 3


def test_p_refinement_study_rates_match_theory():
    results = p_refinement_study(mesh_sizes=(4, 8, 16), orders=(1, 2, 3))
    # Expected L2 rates: P1≈2, P2≈3, P3≈4
    # (docs/science/benchmark_references.md Section 3)
    assert results[1][-1]['l2_rate'] == pytest.approx(2.0, abs=0.15)
    assert results[2][-1]['l2_rate'] == pytest.approx(3.0, abs=0.1)
    assert results[3][-1]['l2_rate'] == pytest.approx(4.0, abs=0.1)


def test_p_refinement_higher_order_is_more_accurate_at_fixed_mesh():
    results = p_refinement_study(mesh_sizes=(8,), orders=(1, 2, 3))
    e1 = results[1][0]['l2_error']
    e2 = results[2][0]['l2_error']
    e3 = results[3][0]['l2_error']
    assert e1 > e2 > e3


# ---------------------------------------------------------------------------
# Convection-diffusion convergence: Galerkin vs. SUPG (Task 4.4)
# ---------------------------------------------------------------------------

def test_convection_diffusion_convergence_structure():
    results = convection_diffusion_convergence(pe_values=(0.1, 100), mesh_sizes=(8, 16))
    assert set(results.keys()) == {0.1, 100}
    for pe, rows in results.items():
        assert len(rows) == 2
        for r in rows:
            assert set(r.keys()) == {
                'n', 'h', 'pe_h',
                'l2_galerkin', 'l2_supg',
                'min_galerkin', 'max_galerkin',
                'min_supg', 'max_supg',
            }


def test_convection_diffusion_low_peclet_both_methods_accurate():
    """At Pe=0.1 (diffusion-dominated), Galerkin and SUPG should agree
    closely and both stay within the exact solution's [0,1] bounds."""
    results = convection_diffusion_convergence(pe_values=(0.1,), mesh_sizes=(8, 32))
    for r in results[0.1]:
        assert r['l2_galerkin'] == pytest.approx(r['l2_supg'], rel=0.05)
        assert r['min_galerkin'] >= -1e-6
        assert r['min_supg'] >= -1e-6


def test_convection_diffusion_high_peclet_supg_is_bounded():
    """At Pe=100 with Pe_h > 1 (coarse mesh), standard Galerkin produces
    spurious undershoot below the exact solution's lower bound (0), while
    SUPG remains bounded in [0,1] -- the key qualitative advantage of SUPG
    stabilization for convection-dominated problems."""
    results = convection_diffusion_convergence(pe_values=(100,), mesh_sizes=(8, 64))
    coarse = results[100][0]
    fine = results[100][1]

    assert coarse['pe_h'] > 1  # convection-dominated on the coarse mesh
    assert coarse['min_galerkin'] < -0.05  # Galerkin undershoots
    assert coarse['min_supg'] > -1e-6  # SUPG stays bounded

    assert fine['pe_h'] < 1  # diffusion-dominated on the fine mesh
    assert fine['min_galerkin'] > -1e-6  # Galerkin also bounded once Pe_h < 1


# ---------------------------------------------------------------------------
# Plotting helpers (DoD: "all convergence plots use log-log scale with
# slope annotation")
# ---------------------------------------------------------------------------

def test_generate_h_refinement_plot_creates_png(tmp_path):
    rows_1d = poisson_1d_convergence(mesh_sizes=(8, 16))
    rows_2d = poisson_2d_convergence(mesh_sizes=(8, 16))
    filename = generate_h_refinement_plot(rows_1d, rows_2d, output_dir=str(tmp_path))
    assert filename == 'h_refinement.png'
    assert (tmp_path / filename).exists()


def test_generate_p_refinement_plot_creates_png(tmp_path):
    p_results = p_refinement_study(mesh_sizes=(4, 8))
    filename = generate_p_refinement_plot(p_results, output_dir=str(tmp_path))
    assert filename == 'p_refinement.png'
    assert (tmp_path / filename).exists()


def test_generate_convection_diffusion_plot_creates_png(tmp_path):
    cd_results = convection_diffusion_convergence(mesh_sizes=(8, 16))
    filename = generate_convection_diffusion_plot(cd_results, output_dir=str(tmp_path))
    assert filename == 'cd_convergence.png'
    assert (tmp_path / filename).exists()


# ---------------------------------------------------------------------------
# De Vahl Davis benchmark (Task 4.5, Goal 7)
# ---------------------------------------------------------------------------

def test_de_vahl_davis_pure_conduction_nusselt_number():
    """For Ra=0 (pure conduction) with hot left wall (u=1), cold right wall
    (u=0), and insulated top/bottom, the exact solution u(x,y)=1-x gives
    Nu = 1.0 exactly -- this verifies the Nusselt-number post-processing."""
    u_sol, mesh, nu_avg = de_vahl_davis_conduction(n=16)
    assert nu_avg == pytest.approx(1.0, abs=1e-8)
    assert mesh.geometry().dim() == 2


def test_de_vahl_davis_solution_matches_linear_profile():
    """The pure-conduction solution should match u(x,y) = 1 - x at several
    sample points, independent of y."""
    u_sol, mesh, _ = de_vahl_davis_conduction(n=16)
    for x, y in [(0.0, 0.25), (0.5, 0.5), (0.5, 0.9), (1.0, 0.1)]:
        assert u_sol(x, y) == pytest.approx(1.0 - x, abs=1e-10)


def test_reference_nusselt_table_contains_ra0_baseline():
    ra0 = next(r for r in REFERENCE_NUSSELT if r['ra'] == 0)
    assert ra0['nu'] == pytest.approx(1.0)
    # De Vahl Davis (1983) reference values for Ra in {1e3..1e6}
    ras = {r['ra'] for r in REFERENCE_NUSSELT}
    assert {1e3, 1e4, 1e5, 1e6}.issubset(ras)


def test_generate_de_vahl_davis_plot_creates_png(tmp_path):
    u_sol, mesh, _ = de_vahl_davis_conduction(n=8)
    filename = generate_de_vahl_davis_plot(u_sol, mesh, output_dir=str(tmp_path))
    assert filename == 'de_vahl_davis.png'
    assert (tmp_path / filename).exists()


# ---------------------------------------------------------------------------
# views.py — formatting helpers
# ---------------------------------------------------------------------------

def test_fmt_rate():
    assert _fmt_rate(None) == '—'
    assert _fmt_rate(2.0) == '2.00'
    assert _fmt_rate(1.999) == '2.00'


def test_fmt_rows():
    rows = poisson_1d_convergence(mesh_sizes=(8, 16))
    formatted = _fmt_rows(rows)
    assert formatted[0]['l2_rate'] == '—'
    assert formatted[1]['l2_rate'] != '—'
    assert 'e' in formatted[0]['l2_error']  # scientific notation
    assert formatted[0]['h'] == '0.1250'


def test_fmt_cd_rows():
    cd_results = convection_diffusion_convergence(pe_values=(100,), mesh_sizes=(8,))
    formatted = _fmt_cd_rows(cd_results[100])
    assert 'e' in formatted[0]['l2_galerkin']
    assert formatted[0]['pe_h'] == '6.250'


def test_check_helper():
    c = _check('name', 'expected', 'observed', True)
    assert c == {'name': 'name', 'expected': 'expected', 'observed': 'observed', 'passed': True}


# ---------------------------------------------------------------------------
# Django views
# ---------------------------------------------------------------------------

def test_validation_url_resolves_to_dashboard(client):
    response = client.get('/Validation/')
    assert response.status_code == 200
    assert 'checks' in response.context


def test_old_steady_state_validation_url_removed(client):
    """The old SteadyStateThermal/validation/ 1D-only check has been
    superseded by the comprehensive Validation dashboard."""
    response = client.get('/SteadyStateThermal/validation/')
    assert response.status_code == 404


def test_dashboard_view_all_checks_pass(client):
    """The dashboard is the project's 'trust anchor' -- every verification
    and benchmark check should pass against its documented expected value."""
    response = client.get('/Validation/')
    assert response.status_code == 200
    checks = response.context['checks']
    assert len(checks) == 8
    for c in checks:
        assert c['passed'], f"check failed: {c['name']} (expected {c['expected']}, got {c['observed']})"


def test_dashboard_check_names_cover_all_dod_items(client):
    response = client.get('/Validation/')
    names = ' '.join(c['name'] for c in response.context['checks'])
    assert '1D Poisson' in names
    assert '2D Poisson' in names
    assert 'p-refinement' in names
    assert 'P2' in names and 'P3' in names
    assert 'SUPG' in names
    assert 'De Vahl Davis' in names


def test_convergence_view(client):
    response = client.get('/Validation/convergence/')
    assert response.status_code == 200
    assert response.context['plot_h'] == 'h_refinement.png'
    assert response.context['plot_p'] == 'p_refinement.png'
    assert response.context['plot_cd'] == 'cd_convergence.png'

    rows_1d = response.context['rows_1d']
    assert len(rows_1d) == 5
    assert rows_1d[0]['l2_rate'] == '—'
    assert rows_1d[-1]['l2_rate'] != '—'

    p_results = dict(response.context['p_results'])
    assert set(p_results.keys()) == {1, 2, 3}

    cd_results = dict(response.context['cd_results'])
    assert set(cd_results.keys()) == {0.1, 100}


def test_benchmark_view(client):
    response = client.get('/Validation/benchmark/')
    assert response.status_code == 200
    assert response.context['plot'] == 'de_vahl_davis.png'
    assert response.context['nu_avg'] == pytest.approx(1.0, abs=1e-6)
    assert response.context['reference'] == REFERENCE_NUSSELT


@pytest.mark.django_db
@pytest.mark.parametrize('path', [
    '/SteadyStateThermal/',
    '/ConvectionDiffusion/',
    '/ConvectionDiffusion/peclet/',
    '/TransientThermal/',
    '/TransientThermal/comparison/',
])
def test_navbar_validation_link_still_resolves_on_other_pages(client, path):
    """Every page's navbar links to {% url 'validation' %}; this must
    continue to resolve (now to the new Validation dashboard) without
    breaking template rendering on any existing page."""
    response = client.get(path)
    assert response.status_code == 200
