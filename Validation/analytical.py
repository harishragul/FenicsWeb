"""Exact (manufactured) solutions used by the Validation module.

Each function returns a FEniCS `Expression` evaluated symbolically (as C++
code, JIT-compiled by FEniCS) so it can be used directly with `errornorm`.
"""
from fenics import Expression


def exact_1d_poisson(degree=5):
    """Exact solution of -u'' = 2 on [0,1], u(0) = u(1) = 0.

    u(x) = x(1 - x)
    """
    return Expression("x[0]*(1 - x[0])", degree=degree)


def exact_2d_poisson(degree=6):
    """Manufactured solution u(x,y) = sin(pi x) sin(pi y) on the unit square,
    u = 0 on the boundary."""
    return Expression("sin(pi*x[0])*sin(pi*x[1])", degree=degree)


def source_2d_poisson(k=1.0, degree=6):
    """Source term f such that -k*lap(u) = f is satisfied by
    u(x,y) = sin(pi x) sin(pi y):

    f(x,y) = 2*pi^2*k*sin(pi x)*sin(pi y)
    """
    return Expression("2*pi*pi*k*sin(pi*x[0])*sin(pi*x[1])", k=k, degree=degree)
