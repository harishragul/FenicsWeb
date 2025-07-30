from fenics import *
def solve_conduction(dimension, mesh, bc_data, f, k):
    """
    Solve heat conduction with support for Dirichlet and convection boundary conditions.
    
    bc_data: dictionary containing boundary condition information
    {
        'left': {'type': 'dirichlet'/'convection', 'value': temp, 'h': coeff, 't_amb': ambient_temp},
        'right': {...},
        'top': {...},
        'bottom': {...},
        'front': {...},
        'back': {...}
    }
    """
    V = FunctionSpace(mesh, "P", 1)
    u = TrialFunction(V)
    v = TestFunction(V)

    # Variational form - conduction term
    a = k * dot(grad(u), grad(v)) * dx
    L = Constant(f) * v * dx

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

    # Define boundary markers
    boundary_markers = MeshFunction("size_t", mesh, mesh.topology().dim()-1)
    boundary_markers.set_all(0)
    
    # Mark boundaries
    boundaries = {}
    if dimension == '1D':
        boundaries = {
            'left': f"near(x[0], {mesh_x_min})",
            'right': f"near(x[0], {mesh_x_max})"
        }
    elif dimension == '2D':
        boundaries = {
            'left': f"near(x[0], {mesh_x_min})",
            'right': f"near(x[0], {mesh_x_max})",
            'bottom': f"near(x[1], {mesh_y_min})",
            'top': f"near(x[1], {mesh_y_max})"
        }
    elif dimension == '3D':
        boundaries = {
            'left': f"near(x[0], {mesh_x_min})",
            'right': f"near(x[0], {mesh_x_max})",
            'bottom': f"near(x[1], {mesh_y_min})",
            'top': f"near(x[1], {mesh_y_max})",
            'front': f"near(x[2], {mesh_z_min})",
            'back': f"near(x[2], {mesh_z_max})"
        }
    
    # Mark boundaries for integration
    marker_id = 1
    tolerance = 1e-14
    
    for boundary_name, boundary_expr in boundaries.items():
        if boundary_name in bc_data and bc_data[boundary_name]['type'] == 'convection':
            if boundary_name == 'left':
                class LeftBoundary(SubDomain):
                    def inside(self, x, on_boundary):
                        return on_boundary and near(x[0], mesh_x_min, tolerance)
                boundary_marker = LeftBoundary()
            elif boundary_name == 'right':
                class RightBoundary(SubDomain):
                    def inside(self, x, on_boundary):
                        return on_boundary and near(x[0], mesh_x_max, tolerance)
                boundary_marker = RightBoundary()
            elif boundary_name == 'bottom' and coordinates.shape[1] >= 2:
                class BottomBoundary(SubDomain):
                    def inside(self, x, on_boundary):
                        return on_boundary and near(x[1], mesh_y_min, tolerance)
                boundary_marker = BottomBoundary()
            elif boundary_name == 'top' and coordinates.shape[1] >= 2:
                class TopBoundary(SubDomain):
                    def inside(self, x, on_boundary):
                        return on_boundary and near(x[1], mesh_y_max, tolerance)
                boundary_marker = TopBoundary()
            elif boundary_name == 'front' and coordinates.shape[1] >= 3:
                class FrontBoundary(SubDomain):
                    def inside(self, x, on_boundary):
                        return on_boundary and near(x[2], mesh_z_min, tolerance)
                boundary_marker = FrontBoundary()
            elif boundary_name == 'back' and coordinates.shape[1] >= 3:
                class BackBoundary(SubDomain):
                    def inside(self, x, on_boundary):
                        return on_boundary and near(x[2], mesh_z_max, tolerance)
                boundary_marker = BackBoundary()
            else:
                continue
                
            boundary_marker.mark(boundary_markers, marker_id)
            bc_data[boundary_name]['marker_id'] = marker_id
            marker_id += 1

    # Define measure for boundary integration
    ds = Measure('ds', domain=mesh, subdomain_data=boundary_markers)
    
    # Boundary conditions
    bcs = []
    
    # Process each boundary
    for boundary_name, boundary_expr in boundaries.items():
        if boundary_name in bc_data:
            bc_info = bc_data[boundary_name]
            
            if bc_info['type'] == 'dirichlet':
                # Dirichlet boundary condition
                bc_value = bc_info['value']
                if bc_value is not None:
                    bcs.append(DirichletBC(V, Constant(bc_value), boundary_expr))
            
            elif bc_info['type'] == 'convection':
                # Convection boundary condition - add to weak form
                h = bc_info['h']
                t_amb = bc_info['t_amb']
                if h is not None and t_amb is not None:
                    marker_id = bc_info['marker_id']
                    # Add convection terms to bilinear and linear forms
                    a += h * u * v * ds(marker_id)
                    L += h * Constant(t_amb) * v * ds(marker_id)

    # Solve
    u_sol = Function(V)
    solve(a == L, u_sol, bcs)

    print(u_sol)

    # Export solution
    solution = u_sol.compute_vertex_values(mesh)
    coordinates = mesh.coordinates()

    solution_list = [[list(coordinates[i]), solution[i]] for i in range(len(solution))]

    return u_sol, solution_list

import os
import matplotlib.pyplot as plt
from FenicsWeb.settings import BASE_DIR
import matplotlib
matplotlib.use("Agg")

STATIC_DIR = os.path.join(BASE_DIR, 'static') 

def generate_solution_plot(solution):
    """Generate a plot of the solution with proper handling for different dimensions"""
    import numpy as np
    
    plt.figure(figsize=(10, 6))
    
    # Get mesh and solution data
    mesh = solution.function_space().mesh()
    coordinates = mesh.coordinates()
    solution_values = solution.compute_vertex_values(mesh)
    
    # Handle different dimensions
    if coordinates.shape[1] == 1:  # 1D problem
        # Sort by x-coordinate for proper line plot
        sorted_indices = np.argsort(coordinates[:, 0])
        x_coords = coordinates[sorted_indices, 0]
        y_values = solution_values[sorted_indices]
        
        plt.plot(x_coords, y_values, 'b-', linewidth=2, marker='o', markersize=4)
        plt.xlabel('X-coordinate')
        plt.ylabel('Temperature (°C)')
        plt.title('Temperature Distribution (1D)')
        plt.grid(True, alpha=0.3)
        
    elif coordinates.shape[1] == 2:  # 2D problem
        try:
            # Try using FEniCS plot function first
            c = plot(solution, cmap='viridis')
            plt.colorbar(c, label='Temperature (°C)')
            plt.xlabel('X-coordinate')
            plt.ylabel('Y-coordinate')
            plt.title('Temperature Distribution (2D)')
        except:
            # Fallback to scatter plot if FEniCS plot fails
            x_coords = coordinates[:, 0]
            y_coords = coordinates[:, 1]
            
            scatter = plt.scatter(x_coords, y_coords, c=solution_values, 
                                cmap='viridis', s=50, edgecolors='black', linewidth=0.5)
            plt.colorbar(scatter, label='Temperature (°C)')
            plt.xlabel('X-coordinate')
            plt.ylabel('Y-coordinate')
            plt.title('Temperature Distribution (2D)')
            
    elif coordinates.shape[1] == 3:  # 3D problem
        # For 3D, create a 2D projection or slice
        x_coords = coordinates[:, 0]
        y_coords = coordinates[:, 1]
        z_coords = coordinates[:, 2]
        
        # Create a scatter plot of the z=0 slice or middle slice
        z_middle = (z_coords.min() + z_coords.max()) / 2
        tolerance = (z_coords.max() - z_coords.min()) * 0.1
        
        mask = np.abs(z_coords - z_middle) < tolerance
        if np.sum(mask) > 0:
            scatter = plt.scatter(x_coords[mask], y_coords[mask], c=solution_values[mask], 
                                cmap='viridis', s=50, edgecolors='black', linewidth=0.5)
            plt.colorbar(scatter, label='Temperature (°C)')
            plt.xlabel('X-coordinate')
            plt.ylabel('Y-coordinate')
            plt.title(f'Temperature Distribution (3D slice at z≈{z_middle:.3f})')
        else:
            # Fallback: plot all points
            scatter = plt.scatter(x_coords, y_coords, c=solution_values, 
                                cmap='viridis', s=30, alpha=0.6)
            plt.colorbar(scatter, label='Temperature (°C)')
            plt.xlabel('X-coordinate')
            plt.ylabel('Y-coordinate')
            plt.title('Temperature Distribution (3D - all points)')
    
    plt.tight_layout()
    
    if not os.path.exists(STATIC_DIR):
        os.makedirs(STATIC_DIR)

    plot_filename = os.path.join(STATIC_DIR, 'solution_plot.png')
    plt.savefig(plot_filename, dpi=150, bbox_inches='tight')
    plt.close()

    return 'solution_plot.png'