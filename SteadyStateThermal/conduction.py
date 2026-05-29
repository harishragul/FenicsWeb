import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from fenics import (
    FunctionSpace, TrialFunction, TestFunction, Function,
    Constant, DirichletBC, MeshFunction, Measure,
    dot, grad, dx,
    solve, near, facets,
)
from FenicsWeb.settings import BASE_DIR

STATIC_DIR = os.path.join(BASE_DIR, 'static')


def _mark_facets(mesh, coords):
    """Return a MeshFunction with facets labelled 1–6 (left/right/bottom/top/front/back)."""
    dim = mesh.geometry().dim()
    markers = MeshFunction("size_t", mesh, dim - 1, 0)

    x_min, x_max = coords[:, 0].min(), coords[:, 0].max()

    if dim == 1:
        for facet in facets(mesh):
            mp = facet.midpoint()
            if near(mp[0], x_min):
                markers[facet] = 1   # left
            elif near(mp[0], x_max):
                markers[facet] = 2   # right

    elif dim == 2:
        y_min, y_max = coords[:, 1].min(), coords[:, 1].max()
        for facet in facets(mesh):
            mp = facet.midpoint()
            if near(mp[0], x_min):
                markers[facet] = 1   # left
            elif near(mp[0], x_max):
                markers[facet] = 2   # right
            elif near(mp[1], y_min):
                markers[facet] = 3   # bottom
            elif near(mp[1], y_max):
                markers[facet] = 4   # top

    elif dim == 3:
        y_min, y_max = coords[:, 1].min(), coords[:, 1].max()
        z_min, z_max = coords[:, 2].min(), coords[:, 2].max()
        for facet in facets(mesh):
            mp = facet.midpoint()
            if near(mp[0], x_min):
                markers[facet] = 1   # left
            elif near(mp[0], x_max):
                markers[facet] = 2   # right
            elif near(mp[1], y_min):
                markers[facet] = 3   # bottom
            elif near(mp[1], y_max):
                markers[facet] = 4   # top
            elif near(mp[2], z_min):
                markers[facet] = 5   # front
            elif near(mp[2], z_max):
                markers[facet] = 6   # back

    return markers


def solve_conduction(dimension, mesh, left_bc, right_bc, top_bc, bottom_bc,
                     front_bc, back_bc, f, k,
                     bc_types=None, neumann_fluxes=None,
                     robin_hs=None, robin_u_infs=None):
    """Solve steady-state heat conduction via FEM; return (u_sol, solution_list).

    Supports Dirichlet, Neumann, and Robin boundary conditions per face.
    bc_types: dict mapping face name to 'dirichlet' | 'neumann' | 'robin'.
              Defaults to all-Dirichlet when None.
    neumann_fluxes: dict of face -> flux value q (used when bc_type=='neumann').
    robin_hs:       dict of face -> convection coefficient h (bc_type=='robin').
    robin_u_infs:   dict of face -> ambient temperature u_inf (bc_type=='robin').
    """
    if bc_types is None:
        bc_types = {}
    if neumann_fluxes is None:
        neumann_fluxes = {}
    if robin_hs is None:
        robin_hs = {}
    if robin_u_infs is None:
        robin_u_infs = {}

    # Default every face to Dirichlet
    faces_1d = ['left', 'right']
    faces_2d = ['left', 'right', 'bottom', 'top']
    faces_3d = ['left', 'right', 'bottom', 'top', 'front', 'back']
    face_map = {
        'left': 1, 'right': 2, 'bottom': 3, 'top': 4, 'front': 5, 'back': 6,
    }
    dirichlet_values = {
        'left': left_bc, 'right': right_bc,
        'top': top_bc, 'bottom': bottom_bc,
        'front': front_bc, 'back': back_bc,
    }

    if dimension == '1D':
        active_faces = faces_1d
    elif dimension == '2D':
        active_faces = faces_2d
    else:
        active_faces = faces_3d

    coordinates = mesh.coordinates()
    markers = _mark_facets(mesh, coordinates)
    ds = Measure('ds', domain=mesh, subdomain_data=markers)

    V = FunctionSpace(mesh, "P", 1)
    u = TrialFunction(V)
    v = TestFunction(V)

    a = k * dot(grad(u), grad(v)) * dx
    L = Constant(f) * v * dx

    bcs = []
    for face in active_faces:
        ftype = bc_types.get(face, 'dirichlet')
        tag = face_map[face]

        if ftype == 'dirichlet':
            x_min = coordinates[:, 0].min()
            x_max = coordinates[:, 0].max()
            # Build subdomain string for DirichletBC
            if face == 'left':
                bc_expr = f"near(x[0], {x_min})"
            elif face == 'right':
                bc_expr = f"near(x[0], {x_max})"
            elif face == 'bottom':
                bc_expr = f"near(x[1], {coordinates[:, 1].min()})"
            elif face == 'top':
                bc_expr = f"near(x[1], {coordinates[:, 1].max()})"
            elif face == 'front':
                bc_expr = f"near(x[2], {coordinates[:, 2].min()})"
            elif face == 'back':
                bc_expr = f"near(x[2], {coordinates[:, 2].max()})"
            bcs.append(DirichletBC(V, Constant(dirichlet_values[face]), bc_expr))

        elif ftype == 'neumann':
            # -k ∂u/∂n = q  →  add q*v*ds to L (natural BC, Neumann flux)
            q = neumann_fluxes.get(face, 0.0)
            L = L + Constant(q) * v * ds(tag)

        elif ftype == 'robin':
            # -k ∂u/∂n = h(u - u_inf)  →  h*u*v*ds added to a, h*u_inf*v*ds added to L
            h = robin_hs.get(face, 0.0)
            u_inf = robin_u_infs.get(face, 0.0)
            a = a + Constant(h) * u * v * ds(tag)
            L = L + Constant(h) * Constant(u_inf) * v * ds(tag)

    u_sol = Function(V)
    solve(a == L, u_sol, bcs)

    solution = u_sol.compute_vertex_values(mesh)
    coordinates = mesh.coordinates()
    solution_list = [[list(coordinates[i]), float(solution[i])] for i in range(len(solution))]

    return u_sol, solution_list


def generate_solution_plot(u_sol, mesh, output_dir=None):
    """Plot the temperature field; dispatch on mesh dimension. Return filename."""
    dim = mesh.geometry().dim()
    target_dir = output_dir if output_dir is not None else STATIC_DIR
    os.makedirs(target_dir, exist_ok=True)
    plot_filename = os.path.join(target_dir, 'solution_plot.png')

    plt.figure(figsize=(8, 6))

    if dim == 1:
        coords = mesh.coordinates().flatten()
        vals = u_sol.compute_vertex_values(mesh)
        order = coords.argsort()
        plt.plot(coords[order], vals[order], 'b-', linewidth=2)
        plt.xlabel('x')
        plt.ylabel('Temperature u(x)')
        plt.title('Temperature Distribution (1D)')
        plt.grid(True)

    elif dim == 2:
        from fenics import plot as fenics_plot
        p = fenics_plot(u_sol)
        plt.colorbar(p, label='Temperature (K)')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('Temperature Distribution (2D)')

    elif dim == 3:
        # Render three cross-section slices via mesh vertex projection
        import numpy as np
        coords = mesh.coordinates()
        vals = u_sol.compute_vertex_values(mesh)
        x, y, z = coords[:, 0], coords[:, 1], coords[:, 2]
        x_mid = (x.max() + x.min()) / 2
        y_mid = (y.max() + y.min()) / 2
        z_mid = (z.max() + z.min()) / 2
        tol = (x.max() - x.min()) / mesh.num_cells() ** (1 / 3)

        fig, axes = plt.subplots(1, 3, figsize=(15, 5))
        slices = [
            ('z = {:.2f}'.format(z_mid), x[abs(z - z_mid) < tol], y[abs(z - z_mid) < tol], vals[abs(z - z_mid) < tol], 'x', 'y'),
            ('y = {:.2f}'.format(y_mid), x[abs(y - y_mid) < tol], z[abs(y - y_mid) < tol], vals[abs(y - y_mid) < tol], 'x', 'z'),
            ('x = {:.2f}'.format(x_mid), y[abs(x - x_mid) < tol], z[abs(x - x_mid) < tol], vals[abs(x - x_mid) < tol], 'y', 'z'),
        ]
        for ax, (title, px, py, pv, xl, yl) in zip(axes, slices):
            if len(px) > 2:
                sc = ax.scatter(px, py, c=pv, cmap='viridis', s=20)
                plt.colorbar(sc, ax=ax, label='T (K)')
            ax.set_title(title)
            ax.set_xlabel(xl)
            ax.set_ylabel(yl)
        plt.suptitle('Temperature Distribution (3D cross-sections)')
        plt.tight_layout()
        plt.savefig(plot_filename, dpi=100, bbox_inches='tight')
        plt.close()
        return 'solution_plot.png'

    plt.tight_layout()
    plt.savefig(plot_filename, dpi=100, bbox_inches='tight')
    plt.close()
    return 'solution_plot.png'
