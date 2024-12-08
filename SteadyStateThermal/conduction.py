from fenics import *
def solve_conduction(dimension, mesh, left_bc, right_bc, top_bc, bottom_bc, front_bc, back_bc, f, k):
    V = FunctionSpace(mesh, "P", 1)
    u = TrialFunction(V)
    v = TestFunction(V)

    # Variational form
    a = k * dot(grad(u), grad(v)) * dx
    L = Constant(f) * v * dx

    # Boundary conditions
    bcs = []
   # Get mesh coordinates and compute min/max for each direction
    coordinates = mesh.coordinates()

    # For 1D
    if coordinates.shape[1] == 1:
        mesh_x_min, mesh_x_max = coordinates.min(), coordinates.max()

    # For 2D
    elif coordinates.shape[1] == 2:
        mesh_x_min, mesh_x_max = coordinates[:, 0].min(), coordinates[:, 0].max()
        mesh_y_min, mesh_y_max = coordinates[:, 1].min(), coordinates[:, 1].max()

    # For 3D
    elif coordinates.shape[1] == 3:
        mesh_x_min, mesh_x_max = coordinates[:, 0].min(), coordinates[:, 0].max()
        mesh_y_min, mesh_y_max = coordinates[:, 1].min(), coordinates[:, 1].max()
        mesh_z_min, mesh_z_max = coordinates[:, 2].min(), coordinates[:, 2].max()
        
    # Adjust boundary conditions based on the actual mesh domain
    if dimension == '1D':
        bcs.append(DirichletBC(V, Constant(left_bc), f"near(x[0], {mesh_x_min})"))
        bcs.append(DirichletBC(V, Constant(right_bc), f"near(x[0], {mesh_x_max})"))
    elif dimension == '2D':
        bcs.append(DirichletBC(V, Constant(left_bc), f"near(x[0], {mesh_x_min})"))
        bcs.append(DirichletBC(V, Constant(right_bc), f"near(x[0], {mesh_x_max})"))
        bcs.append(DirichletBC(V, Constant(bottom_bc), f"near(x[1], {mesh_y_min})"))
        bcs.append(DirichletBC(V, Constant(top_bc), f"near(x[1], {mesh_y_max})"))
    elif dimension == '3D':
        bcs.append(DirichletBC(V, Constant(left_bc), f"near(x[0], {mesh_x_min})"))
        bcs.append(DirichletBC(V, Constant(right_bc), f"near(x[0], {mesh_x_max})"))
        bcs.append(DirichletBC(V, Constant(bottom_bc), f"near(x[1], {mesh_y_min})"))
        bcs.append(DirichletBC(V, Constant(top_bc), f"near(x[1], {mesh_y_max})"))
        bcs.append(DirichletBC(V, Constant(front_bc), f"near(x[2], {mesh_z_min})"))
        bcs.append(DirichletBC(V, Constant(back_bc), f"near(x[2], {mesh_z_max})"))
        

    # Solve
    u_sol = Function(V)
    solve(a == L, u_sol, bcs)

    # Export solution
    solution = u_sol.vector().get_local()
    coordinates = mesh.coordinates()

    return u_sol, solution, coordinates

import os
import matplotlib.pyplot as plt
from FenicsWeb.settings import BASE_DIR
import matplotlib
matplotlib.use("Agg")

STATIC_DIR = os.path.join(BASE_DIR, 'static') 

def generate_solution_plot(mesh):
    plot(mesh)

    if not os.path.exists(STATIC_DIR):
        os.makedirs(STATIC_DIR)

    plot_filename = os.path.join(STATIC_DIR, 'solution_plot.png')
    plt.savefig(plot_filename)
    plt.close()

    return 'solution_plot.png'