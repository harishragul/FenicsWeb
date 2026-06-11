import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from fenics import (
    FunctionSpace, TrialFunction, TestFunction, Function,
    Constant, DirichletBC, MeshFunction, CellDiameter,
    dot, grad, dx, sqrt,
    solve, facets, errornorm, Expression,
    UnitIntervalMesh, DOLFIN_EPS,
)
from ufl import tanh, conditional, gt
from FenicsWeb.settings import BASE_DIR

STATIC_DIR = os.path.join(BASE_DIR, 'static')


def compute_peclet_number(mesh, k, bx, by=0.0):
    """Return the (global) cell Peclet number Pe_h = |b| * h / (2k).

    Uses the largest cell diameter in the mesh as the representative
    length scale h. Pe_h < 1 indicates a diffusion-dominated regime where
    standard Galerkin is stable; Pe_h > 1 indicates a convection-dominated
    regime where standard Galerkin produces spurious oscillations and SUPG
    stabilization is required.
    """
    if k <= 0:
        return float('inf')
    b_norm = (bx ** 2 + by ** 2) ** 0.5
    h = mesh.hmax()
    return b_norm * h / (2.0 * k)


def solve_convection_diffusion(dimension, mesh, k, bx, by, f,
                                left_bc=0.0, right_bc=1.0,
                                bottom_bc=0.0, top_bc=0.0,
                                method='galerkin'):
    """Solve the steady convection-diffusion equation -k*grad^2(u) + b.grad(u) = f.

    Dirichlet boundary conditions are applied on every boundary of the mesh
    (left/right for 1D; left/right/bottom/top for 2D), using the mesh's own
    coordinate bounds.

    method: 'galerkin' for the standard Galerkin form, or 'supg' for the
    SUPG-stabilized form (residual-weighted streamline test function).

    Returns (u_sol, solution_list) where solution_list mirrors the format
    used by SteadyStateThermal.solve_conduction.
    """
    V = FunctionSpace(mesh, "P", 1)
    u = TrialFunction(V)
    v = TestFunction(V)

    b = Constant((bx,)) if dimension == '1D' else Constant((bx, by))

    coords = mesh.coordinates()
    x_min, x_max = float(coords[:, 0].min()), float(coords[:, 0].max())

    bcs = [
        DirichletBC(V, Constant(left_bc), f"near(x[0], {x_min})"),
        DirichletBC(V, Constant(right_bc), f"near(x[0], {x_max})"),
    ]
    if dimension == '2D':
        y_min, y_max = float(coords[:, 1].min()), float(coords[:, 1].max())
        bcs.append(DirichletBC(V, Constant(bottom_bc), f"near(x[1], {y_min})"))
        bcs.append(DirichletBC(V, Constant(top_bc), f"near(x[1], {y_max})"))

    if method == 'supg':
        h = CellDiameter(mesh)
        b_norm = sqrt(dot(b, b))
        Pe = b_norm * h / (2 * k)
        # tau = h/(2|b|) * (coth(Pe) - 1/Pe), with the small-Pe (or b=0)
        # limit tau -> h^2/(12k) handled via `conditional` to avoid 0/0.
        tau = conditional(
            gt(Pe, 1e-8),
            h / (2 * b_norm + DOLFIN_EPS) * (1 / tanh(Pe) - 1 / Pe),
            h * h / (12 * k),
        )
        v_test = v + tau * dot(b, grad(v))
    else:
        v_test = v

    a = k * dot(grad(u), grad(v_test)) * dx + dot(b, grad(u)) * v_test * dx
    L = Constant(f) * v_test * dx

    u_sol = Function(V)
    solve(a == L, u_sol, bcs)

    values = u_sol.compute_vertex_values(mesh)
    solution_list = [[list(coords[i]), float(values[i])] for i in range(len(values))]

    return u_sol, solution_list


def mark_inflow_outflow(mesh, bx, by=0.0):
    """Mark exterior facets according to the sign of b.n (the velocity field
    dotted with the outward normal).

    Returns a MeshFunction over facets with values:
        1 = inflow         (b.n < 0) — Dirichlet (prescribed) boundary
        2 = outflow        (b.n > 0) — natural (zero-flux) boundary
        3 = characteristic (b.n ~ 0) — no-flux wall
        0 = interior facet
    """
    dim = mesh.topology().dim()
    mesh.init(dim - 1, dim)
    markers = MeshFunction("size_t", mesh, dim - 1, 0)
    b_vec = np.array([bx, by]) if dim == 2 else np.array([bx])

    for facet in facets(mesh):
        if not facet.exterior():
            continue
        n = facet.normal()
        n_arr = np.array([n.x(), n.y()]) if dim == 2 else np.array([n.x()])
        bn = float(np.dot(b_vec, n_arr))
        if bn < -DOLFIN_EPS:
            markers[facet] = 1   # inflow
        elif bn > DOLFIN_EPS:
            markers[facet] = 2   # outflow
        else:
            markers[facet] = 3   # characteristic / no-flux

    return markers


def exact_1d_solution(Pe, degree=5):
    """Exact solution of -u'' + Pe*u' = 0 on [0,1], u(0)=0, u(1)=1.

    u(x) = (exp(Pe*x) - 1) / (exp(Pe) - 1). For Pe -> 0 this degenerates to
    the linear profile u(x) = x.
    """
    if abs(Pe) < 1e-8:
        return Expression("x[0]", degree=degree)
    return Expression("(exp(Pe*x[0]) - 1) / (exp(Pe) - 1)", Pe=Pe, degree=degree)


def peclet_study_1d(pe_values, n=20):
    """Run the 1D benchmark -u'' + Pe*u' = 0, u(0)=0, u(1)=1 for each value of
    Pe in `pe_values`, with both standard Galerkin and SUPG, on a fixed mesh
    of `n` cells (k = 1/Pe, b = 1).

    Returns a list of dicts (one per Pe) with the nodal solutions, the exact
    solution, and the L2 error for both methods.
    """
    mesh = UnitIntervalMesh(n)
    coords = mesh.coordinates().flatten()
    order = coords.argsort()
    x_sorted = coords[order]

    results = []
    for Pe in pe_values:
        k = 1.0 / Pe
        bx = 1.0

        u_g, _ = solve_convection_diffusion('1D', mesh, k, bx, 0.0, 0.0,
                                             left_bc=0.0, right_bc=1.0,
                                             method='galerkin')
        u_s, _ = solve_convection_diffusion('1D', mesh, k, bx, 0.0, 0.0,
                                             left_bc=0.0, right_bc=1.0,
                                             method='supg')
        u_exact = exact_1d_solution(Pe)

        vals_g = u_g.compute_vertex_values(mesh)[order]
        vals_s = u_s.compute_vertex_values(mesh)[order]
        vals_e = np.array([u_exact(xi) for xi in x_sorted])

        results.append({
            'pe': Pe,
            'cell_pe': bx * mesh.hmax() / (2 * k),
            'x': x_sorted.tolist(),
            'u_exact': vals_e.tolist(),
            'u_galerkin': vals_g.tolist(),
            'u_supg': vals_s.tolist(),
            'l2_galerkin': errornorm(u_exact, u_g, 'L2'),
            'l2_supg': errornorm(u_exact, u_s, 'L2'),
        })

    return results


def generate_convection_solution_plot(u_sol, mesh, bx, by, label='', output_dir=None):
    """Plot the convection-diffusion solution field. For 2D meshes the
    velocity field b=(bx,by) is overlaid as a quiver plot. Returns the
    output filename (relative to output_dir / STATIC_DIR).
    """
    dim = mesh.geometry().dim()
    target_dir = output_dir if output_dir is not None else STATIC_DIR
    os.makedirs(target_dir, exist_ok=True)
    suffix = f"_{label}" if label else ""
    plot_filename = os.path.join(target_dir, f'convection_solution{suffix}.png')

    plt.figure(figsize=(8, 6))

    if dim == 1:
        coords = mesh.coordinates().flatten()
        vals = u_sol.compute_vertex_values(mesh)
        order = coords.argsort()
        plt.plot(coords[order], vals[order], 'b-', linewidth=2)
        plt.xlabel('x')
        plt.ylabel('u(x)')
        plt.grid(True)
    else:
        from fenics import plot as fenics_plot
        p = fenics_plot(u_sol)
        plt.colorbar(p, label='u')
        if abs(bx) > 0 or abs(by) > 0:
            coords = mesh.coordinates()
            x, y = coords[:, 0], coords[:, 1]
            plt.quiver(x, y, bx * np.ones_like(x), by * np.ones_like(y),
                       color='white', alpha=0.6)
        plt.xlabel('x')
        plt.ylabel('y')

    title = 'Convection-Diffusion Solution'
    if label:
        title += f' ({label.upper()})'
    plt.title(title)

    plt.tight_layout()
    plt.savefig(plot_filename, dpi=100, bbox_inches='tight')
    plt.close()
    return os.path.basename(plot_filename)


def generate_peclet_plot(results, output_dir=None):
    """Generate a side-by-side comparison of the exact, Galerkin, and SUPG
    solutions for each Peclet number in `results` (as returned by
    peclet_study_1d). Returns the output filename.
    """
    target_dir = output_dir if output_dir is not None else STATIC_DIR
    os.makedirs(target_dir, exist_ok=True)
    plot_filename = os.path.join(target_dir, 'peclet_study.png')

    n = len(results)
    fig, axes = plt.subplots(1, n, figsize=(5 * n, 4), sharey=True)
    if n == 1:
        axes = [axes]

    for ax, r in zip(axes, results):
        x = r['x']
        ax.plot(x, r['u_exact'], 'k-', label='Exact', linewidth=2)
        ax.plot(x, r['u_galerkin'], 'r--o', label='Galerkin', markersize=3)
        ax.plot(x, r['u_supg'], 'b-.s', label='SUPG', markersize=3)
        ax.set_title(f"Pe = {r['pe']:.0f}  (Pe_h = {r['cell_pe']:.2f})")
        ax.set_xlabel('x')
        ax.legend(fontsize=8)
        ax.grid(True)
    axes[0].set_ylabel('u(x)')

    plt.suptitle('Galerkin vs. SUPG — 1D Convection-Diffusion Peclet Study')
    plt.tight_layout()
    plt.savefig(plot_filename, dpi=100, bbox_inches='tight')
    plt.close()
    return os.path.basename(plot_filename)
