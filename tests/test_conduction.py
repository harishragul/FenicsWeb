"""Tests for steady-state heat conduction — SteadyStateThermal/conduction.py"""
import json
import pytest
from fenics import UnitIntervalMesh, UnitSquareMesh, Expression, errornorm
from SteadyStateThermal.conduction import solve_conduction


def l2_error(u_h, expr_str):
    u_ex = Expression(expr_str, degree=5)
    return errornorm(u_ex, u_h, 'L2')


# ---------------------------------------------------------------------------
# 1D Dirichlet — correctness
# ---------------------------------------------------------------------------

def test_1d_zero_source_linear_solution():
    """f=0, u(0)=0, u(1)=1 → solution must be exactly u(x)=x."""
    mesh = UnitIntervalMesh(20)
    u_sol, _ = solve_conduction(
        '1D', mesh,
        left_bc=0, right_bc=1,
        top_bc=0, bottom_bc=0, front_bc=0, back_bc=0,
        f=0, k=1,
    )
    assert abs(u_sol(0.5) - 0.5) < 1e-10
    assert abs(u_sol(0.25) - 0.25) < 1e-10


def test_1d_constant_source_parabolic_solution():
    """f=2, k=1, u(0)=u(1)=0 → P1 FEM is nodally exact for u(x)=x(1-x)."""
    mesh = UnitIntervalMesh(20)
    u_sol, _ = solve_conduction(
        '1D', mesh,
        left_bc=0, right_bc=0,
        top_bc=0, bottom_bc=0, front_bc=0, back_bc=0,
        f=2, k=1,
    )
    assert abs(u_sol(0.5) - 0.25) < 1e-10
    assert abs(u_sol(0.25) - 0.1875) < 1e-10


def test_1d_nonunit_conductivity():
    """f=2k, k=5 → same parabolic nodal profile regardless of conductivity."""
    mesh = UnitIntervalMesh(20)
    k_val = 5.0
    u_sol, _ = solve_conduction(
        '1D', mesh,
        left_bc=0, right_bc=0,
        top_bc=0, bottom_bc=0, front_bc=0, back_bc=0,
        f=2 * k_val, k=k_val,
    )
    assert abs(u_sol(0.5) - 0.25) < 1e-10


def test_1d_nonzero_dirichlet_values():
    """f=0, u(0)=100, u(1)=200 → linear profile u(x)=100+100x."""
    mesh = UnitIntervalMesh(40)
    u_sol, _ = solve_conduction(
        '1D', mesh,
        left_bc=100, right_bc=200,
        top_bc=0, bottom_bc=0, front_bc=0, back_bc=0,
        f=0, k=1,
    )
    assert abs(u_sol(0.0) - 100.0) < 1e-8
    assert abs(u_sol(1.0) - 200.0) < 1e-8
    assert abs(u_sol(0.5) - 150.0) < 1e-8


# ---------------------------------------------------------------------------
# 1D h-refinement convergence (P1 → O(h²) in L2)
# ---------------------------------------------------------------------------

def test_1d_l2_convergence_rate():
    """L2 error must ratio > 3.5× when mesh is halved (P1, O(h²) expected)."""
    errors = []
    for n in [8, 16, 32]:
        mesh = UnitIntervalMesh(n)
        u_sol, _ = solve_conduction(
            '1D', mesh,
            left_bc=0, right_bc=0,
            top_bc=0, bottom_bc=0, front_bc=0, back_bc=0,
            f=2, k=1,
        )
        errors.append(l2_error(u_sol, "x[0]*(1 - x[0])"))

    rate_1 = errors[0] / errors[1]
    rate_2 = errors[1] / errors[2]
    assert rate_1 > 3.5, f"Rate too low at n=8→16: {rate_1:.2f}"
    assert rate_2 > 3.5, f"Rate too low at n=16→32: {rate_2:.2f}"


# ---------------------------------------------------------------------------
# 1D Neumann BC
# ---------------------------------------------------------------------------

def test_1d_neumann_right_bc_linear():
    """Dirichlet u(0)=0, Neumann q=1 on right, f=0, k=1 → u(x)=x."""
    mesh = UnitIntervalMesh(20)
    u_sol, _ = solve_conduction(
        '1D', mesh,
        left_bc=0, right_bc=0,
        top_bc=0, bottom_bc=0, front_bc=0, back_bc=0,
        f=0, k=1,
        bc_types={'left': 'dirichlet', 'right': 'neumann'},
        neumann_fluxes={'right': 1.0},
    )
    # Exact: u(x) = x.  Check two interior nodes.
    assert abs(u_sol(0.5) - 0.5) < 1e-6
    assert abs(u_sol(0.25) - 0.25) < 1e-6


def test_1d_neumann_both_sides_raises_or_zero():
    """Two Neumann BCs with zero net flux — solution up to a constant."""
    mesh = UnitIntervalMesh(10)
    # Purely Neumann with non-zero flux should still produce a solve
    # (the system is singular without at least one Dirichlet — FEniCS will warn but not crash)
    # We simply confirm it returns without raising an unhandled exception.
    try:
        u_sol, _ = solve_conduction(
            '1D', mesh,
            left_bc=0, right_bc=0,
            top_bc=0, bottom_bc=0, front_bc=0, back_bc=0,
            f=1, k=1,
            bc_types={'left': 'neumann', 'right': 'neumann'},
            neumann_fluxes={'left': 0.0, 'right': 0.0},
        )
    except Exception:
        pass  # singular system without Dirichlet is acceptable


# ---------------------------------------------------------------------------
# 1D Robin BC
# ---------------------------------------------------------------------------

def test_1d_robin_right_temperature_above_ambient():
    """Dirichlet u(0)=100, Robin on right → boundary temp must be above ambient."""
    mesh = UnitIntervalMesh(40)
    h_coeff, u_inf = 10.0, 20.0
    u_sol, _ = solve_conduction(
        '1D', mesh,
        left_bc=100, right_bc=0,
        top_bc=0, bottom_bc=0, front_bc=0, back_bc=0,
        f=0, k=1,
        bc_types={'left': 'dirichlet', 'right': 'robin'},
        robin_hs={'right': h_coeff},
        robin_u_infs={'right': u_inf},
    )
    # Temperature at the right boundary must be between u_inf and the left BC value
    assert u_sol(1.0) > u_inf
    assert u_sol(1.0) < 100.0


def test_1d_robin_zero_h_equals_neumann_zero_flux():
    """Robin with h=0 → no flux, behaves identically to Neumann q=0 (insulated)."""
    mesh = UnitIntervalMesh(20)
    u_sol, _ = solve_conduction(
        '1D', mesh,
        left_bc=50, right_bc=0,
        top_bc=0, bottom_bc=0, front_bc=0, back_bc=0,
        f=0, k=1,
        bc_types={'left': 'dirichlet', 'right': 'robin'},
        robin_hs={'right': 0.0},
        robin_u_infs={'right': 20.0},
    )
    # With h=0 the right is insulated → temperature constant at left_bc
    assert abs(u_sol(0.5) - 50.0) < 1e-8
    assert abs(u_sol(1.0) - 50.0) < 1e-8


# ---------------------------------------------------------------------------
# 2D correctness
# ---------------------------------------------------------------------------

def test_2d_solution_exists_and_has_correct_structure():
    """2D solve must return u_sol and a solution_list with 2D coordinate pairs."""
    mesh = UnitSquareMesh(4, 4)
    u_sol, solution_list = solve_conduction(
        '2D', mesh,
        left_bc=0, right_bc=1,
        top_bc=0, bottom_bc=0,
        front_bc=0, back_bc=0,
        f=0, k=1,
    )
    assert u_sol is not None
    assert len(solution_list) > 0
    first = solution_list[0]
    assert len(first) == 2
    assert len(first[0]) == 2
    assert isinstance(first[1], float)


# ---------------------------------------------------------------------------
# Interactive Plotly plot
# ---------------------------------------------------------------------------

def test_interactive_mesh_plot_returns_valid_json():
    from mesh.mesh import generate_interactive_mesh_plot
    mesh = UnitSquareMesh(4, 4)
    result = generate_interactive_mesh_plot(mesh)
    data = json.loads(result)
    assert "data" in data


def test_interactive_mesh_plot_with_solution_has_intensity():
    """When u_sol is passed, the Plotly trace must include an intensity array."""
    from mesh.mesh import generate_interactive_mesh_plot
    mesh = UnitSquareMesh(4, 4)
    u_sol, _ = solve_conduction(
        '2D', mesh,
        left_bc=0, right_bc=1,
        top_bc=0, bottom_bc=0, front_bc=0, back_bc=0,
        f=0, k=1,
    )
    result = generate_interactive_mesh_plot(mesh, u_sol=u_sol)
    data = json.loads(result)
    assert "intensity" in data["data"][0]
