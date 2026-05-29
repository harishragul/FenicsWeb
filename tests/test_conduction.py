"""Tests for steady-state heat conduction — SteadyStateThermal/conduction.py"""
import pytest
from fenics import (
    UnitIntervalMesh, UnitSquareMesh,
    Expression, errornorm,
)
from SteadyStateThermal.conduction import solve_conduction


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def l2_error(u_h, u_exact_expr):
    # errornorm computes ||u_h - u_exact||_L2 against the symbolic expression,
    # avoiding the mesh-association pitfall of manually assembling with dx.
    u_exact = Expression(u_exact_expr, degree=5)
    return errornorm(u_exact, u_h, 'L2')


# ---------------------------------------------------------------------------
# 1D correctness tests
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
    """f=2, k=1, u(0)=u(1)=0 → P1 FEM is exact at mesh nodes for u(x)=x(1-x)."""
    mesh = UnitIntervalMesh(20)
    u_sol, _ = solve_conduction(
        '1D', mesh,
        left_bc=0, right_bc=0,
        top_bc=0, bottom_bc=0, front_bc=0, back_bc=0,
        f=2, k=1,
    )
    # P1 Galerkin is exact at nodes for quadratic u=x(1-x) — verify two nodes
    assert abs(u_sol(0.5) - 0.5 * 0.5) < 1e-10      # node at x=0.5: u=0.25
    assert abs(u_sol(0.25) - 0.25 * 0.75) < 1e-10   # node at x=0.25: u=0.1875


def test_1d_nonunit_conductivity():
    """f=2k, k=5 → same parabolic nodal solution regardless of conductivity."""
    mesh = UnitIntervalMesh(20)
    k_val = 5.0
    u_sol, _ = solve_conduction(
        '1D', mesh,
        left_bc=0, right_bc=0,
        top_bc=0, bottom_bc=0, front_bc=0, back_bc=0,
        f=2 * k_val, k=k_val,
    )
    # Scaling k and f together leaves u=x(1-x) unchanged
    assert abs(u_sol(0.5) - 0.25) < 1e-10
    assert abs(u_sol(0.75) - 0.75 * 0.25) < 1e-10


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
# 1D h-refinement convergence test (P1 → O(h²) in L2)
# ---------------------------------------------------------------------------

def test_1d_l2_convergence_rate():
    """L2 error must halve at least 3.5× when mesh is halved (P1, O(h²))."""
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

    rate_1 = errors[0] / errors[1]   # n=8 → n=16, expect ≈ 4
    rate_2 = errors[1] / errors[2]   # n=16 → n=32, expect ≈ 4
    assert rate_1 > 3.5, f"Convergence rate too low: {rate_1:.2f}"
    assert rate_2 > 3.5, f"Convergence rate too low: {rate_2:.2f}"


# ---------------------------------------------------------------------------
# 2D correctness test
# ---------------------------------------------------------------------------

def test_2d_zero_source_solution_exists():
    """2D conduction must return a solution without crashing."""
    mesh = UnitSquareMesh(8, 8)
    u_sol, solution_list = solve_conduction(
        '2D', mesh,
        left_bc=0, right_bc=1,
        top_bc=0, bottom_bc=0,
        front_bc=0, back_bc=0,
        f=0, k=1,
    )
    assert u_sol is not None
    assert len(solution_list) > 0


def test_2d_solution_list_structure():
    """Each entry in solution_list must be [[x, y], temperature_value]."""
    mesh = UnitSquareMesh(4, 4)
    _, solution_list = solve_conduction(
        '2D', mesh,
        left_bc=0, right_bc=1,
        top_bc=0, bottom_bc=0,
        front_bc=0, back_bc=0,
        f=0, k=1,
    )
    first = solution_list[0]
    assert len(first) == 2         # [coords, value]
    assert len(first[0]) == 2      # 2D coords [x, y]
    assert isinstance(first[1], float)
