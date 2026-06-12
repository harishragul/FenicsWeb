"""Benchmark problems for the Validation module (Phase 4 task 4.5,
Goal 7 of the project specification).

Implements the De Vahl Davis differentially-heated cavity in its
pure-conduction (Ra -> 0) limit: hot left wall (u=1), cold right wall (u=0),
insulated (zero-flux / natural Neumann) top and bottom walls. The average
Nusselt number along the hot wall,

    Nu = -∫_{x=0} ∂u/∂x dy

is compared against the De Vahl Davis (1983) reference table. For Ra=0 the
exact solution is u(x,y) = 1 - x (independent of y), so the exact Nusselt
number is Nu = 1.0 -- this serves as a verification target for the Nusselt
number computation itself, since FenicsWeb solves only the steady thermal
diffusion equation (the Ra=0 limit of the full natural-convection benchmark).
"""
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from fenics import (
    Constant, DirichletBC, Function, FunctionSpace, Measure, MeshFunction,
    TestFunction, TrialFunction, UnitSquareMesh, assemble, dot, dx, facets,
    grad, near, plot, solve,
)

from FenicsWeb.settings import BASE_DIR

STATIC_DIR = os.path.join(BASE_DIR, 'static')

# De Vahl Davis (1983) reference average Nusselt numbers for the
# differentially-heated square cavity. FenicsWeb solves only the thermal
# diffusion problem, i.e. the Ra=0 (pure conduction) limit.
REFERENCE_NUSSELT = [
    {'ra': 0, 'nu': 1.000, 'note': 'Pure conduction (this solver)'},
    {'ra': 1e3, 'nu': 1.118, 'note': 'De Vahl Davis (1983)'},
    {'ra': 1e4, 'nu': 2.243, 'note': 'De Vahl Davis (1983)'},
    {'ra': 1e5, 'nu': 4.519, 'note': 'De Vahl Davis (1983)'},
    {'ra': 1e6, 'nu': 8.800, 'note': 'De Vahl Davis (1983)'},
]


def de_vahl_davis_conduction(n=32):
    """Solve the pure-conduction (Ra=0) limit of the De Vahl Davis cavity on
    the unit square: hot left wall u=1, cold right wall u=0, insulated
    (natural Neumann) top/bottom.

    Returns (u_sol, mesh, nu_avg) where nu_avg is the average Nusselt number
    along the hot (left) wall.
    """
    mesh = UnitSquareMesh(n, n)
    V = FunctionSpace(mesh, "P", 1)
    u = TrialFunction(V)
    v = TestFunction(V)

    bcs = [
        DirichletBC(V, Constant(1.0), "near(x[0], 0.0)"),
        DirichletBC(V, Constant(0.0), "near(x[0], 1.0)"),
    ]
    a = dot(grad(u), grad(v)) * dx
    L = Constant(0.0) * v * dx
    u_sol = Function(V)
    solve(a == L, u_sol, bcs)

    markers = MeshFunction("size_t", mesh, mesh.topology().dim() - 1, 0)
    for facet in facets(mesh):
        if near(facet.midpoint()[0], 0.0):
            markers[facet] = 1
    ds = Measure('ds', domain=mesh, subdomain_data=markers)

    nu_avg = -assemble(grad(u_sol)[0] * ds(1))
    return u_sol, mesh, float(nu_avg)


def generate_de_vahl_davis_plot(u_sol, mesh, output_dir=None):
    """Plot the pure-conduction temperature field of the De Vahl Davis
    cavity. Returns the output filename."""
    target_dir = output_dir if output_dir is not None else STATIC_DIR
    os.makedirs(target_dir, exist_ok=True)
    plot_filename = os.path.join(target_dir, 'de_vahl_davis.png')

    fig, ax = plt.subplots(figsize=(6, 4.8))
    c = plot(u_sol)
    plt.colorbar(c, ax=ax, label='u (dimensionless temperature)')
    ax.set_title('De Vahl Davis cavity — pure conduction (Ra = 0)', fontsize=11)
    ax.set_xlabel('x')
    ax.set_ylabel('y')

    plt.tight_layout()
    plt.savefig(plot_filename, dpi=120)
    plt.close(fig)
    return 'de_vahl_davis.png'
