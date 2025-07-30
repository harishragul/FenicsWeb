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
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import numpy as np

STATIC_DIR = os.path.join(BASE_DIR, 'static') 

def generate_solution_plot(solution):
    """Generate an interactive plot of the solution with proper handling for different dimensions"""
    
    # Get mesh and solution data
    mesh = solution.function_space().mesh()
    coordinates = mesh.coordinates()
    solution_values = solution.compute_vertex_values(mesh)
    
    # Handle different dimensions
    if coordinates.shape[1] == 1:  # 1D problem
        sorted_indices = np.argsort(coordinates[:, 0])
        x_coords = coordinates[sorted_indices, 0]
        y_values = solution_values[sorted_indices]
        
        # Calculate heat flow (negative gradient)
        grad_T = np.gradient(y_values, x_coords)
        q = -grad_T  # Heat flux direction (Fourier's law: q = -k * dT/dx)
        
        # Create interactive 1D plot with three subplots
        fig = make_subplots(
            rows=3, cols=1,
            subplot_titles=('Temperature Distribution', 'Temperature Gradient (dT/dx)', 'Heat Flow (-k*dT/dx)'),
            vertical_spacing=0.12,
            row_heights=[0.4, 0.3, 0.3]
        )
        
        # Temperature plot
        fig.add_trace(
            go.Scatter(
                x=x_coords, 
                y=y_values,
                mode='lines+markers',
                name='Temperature',
                line=dict(color='blue', width=3),
                marker=dict(size=6),
                hovertemplate='X: %{x:.3f}<br>Temperature: %{y:.3f}°C<extra></extra>'
            ),
            row=1, col=1
        )
        
        # Temperature gradient plot
        fig.add_trace(
            go.Scatter(
                x=x_coords,
                y=grad_T,
                mode='lines+markers',
                name='dT/dx',
                line=dict(color='green', width=2),
                marker=dict(size=4),
                hovertemplate='X: %{x:.3f}<br>Gradient: %{y:.3f}°C/m<extra></extra>'
            ),
            row=2, col=1
        )
        
        # Heat flow plot
        fig.add_trace(
            go.Scatter(
                x=x_coords,
                y=q,
                mode='lines+markers',
                name='Heat Flow',
                line=dict(color='red', width=2),
                marker=dict(size=4),
                hovertemplate='X: %{x:.3f}<br>Heat Flow: %{y:.3f} W/m²<extra></extra>',
                fill='tonexty' if np.any(q < 0) else 'tozeroy',
                fillcolor='rgba(255,0,0,0.2)'
            ),
            row=3, col=1
        )
        
        # Add arrows to show heat flow direction on temperature plot
        skip = max(1, len(x_coords)//15)
        arrow_scale = (max(y_values) - min(y_values)) * 0.1  # Scale arrows relative to temp range
        
        for i in range(0, len(x_coords), skip):
            if abs(q[i]) > 1e-10:  # Only show arrows where there's significant flow
                # Arrow pointing in direction of heat flow
                arrow_dir = 1 if q[i] > 0 else -1
                x_pos = x_coords[i]
                y_pos = y_values[i]
                
                # Add arrow as annotation
                fig.add_annotation(
                    x=x_pos, y=y_pos,
                    ax=x_pos + arrow_dir * 0.05 * (max(x_coords) - min(x_coords)),
                    ay=y_pos,
                    arrowhead=2, arrowsize=1.5, arrowwidth=3,
                    arrowcolor="red", opacity=0.8,
                    showarrow=True, axref="x", ayref="y",
                    row=1, col=1
                )
        
        fig.update_xaxes(title_text="X-coordinate", row=1, col=1)
        fig.update_xaxes(title_text="X-coordinate", row=2, col=1)
        fig.update_xaxes(title_text="X-coordinate", row=3, col=1)
        fig.update_yaxes(title_text="Temperature (°C)", row=1, col=1)
        fig.update_yaxes(title_text="dT/dx (°C/m)", row=2, col=1)
        fig.update_yaxes(title_text="Heat Flow (W/m²)", row=3, col=1)
        
        # Add horizontal line at zero for gradient and heat flow
        fig.add_hline(y=0, line_dash="dash", line_color="gray", opacity=0.5, row=2, col=1)
        fig.add_hline(y=0, line_dash="dash", line_color="gray", opacity=0.5, row=3, col=1)
        
        fig.update_layout(
            title="Interactive 1D Temperature Distribution and Heat Flow Analysis",
            height=800,
            showlegend=True
        )
        
    elif coordinates.shape[1] == 2:  # 2D problem
        x_coords = coordinates[:, 0]
        y_coords = coordinates[:, 1]
        
        # Calculate heat flow vectors
        V = solution.function_space()
        mesh = V.mesh()
        grad_T = project(grad(solution), VectorFunctionSpace(mesh, 'P', 1))
        grad_vals = grad_T.compute_vertex_values(mesh)
        N = coordinates.shape[0]
        grad_x = -grad_vals[:N]  # negative gradient (heat flow x)
        grad_y = -grad_vals[N:]  # negative gradient (heat flow y)
        
        # Create 2D interactive plot
        fig = go.Figure()
        
        # Temperature contour/scatter plot
        fig.add_trace(
            go.Scatter(
                x=x_coords,
                y=y_coords,
                mode='markers',
                marker=dict(
                    size=8,
                    color=solution_values,
                    colorscale='Viridis',
                    colorbar=dict(title="Temperature (°C)", x=1.15),
                    line=dict(width=1, color='black')
                ),
                name='Temperature',
                hovertemplate='X: %{x:.3f}<br>Y: %{y:.3f}<br>Temperature: %{marker.color:.3f}°C<extra></extra>'
            )
        )
        
        # Add heat flow vectors using proper quiver-like visualization
        skip = max(1, N//25)  # Subsample for better visualization
        scale_factor = 0.03  # Scale arrows appropriately
        
        # Get subsampled coordinates and heat flow vectors
        x_sub = x_coords[::skip]
        y_sub = y_coords[::skip]
        u_sub = grad_x[::skip] * scale_factor
        v_sub = grad_y[::skip] * scale_factor
        
        # Add a single legend entry for heat flow first
        fig.add_trace(
            go.Scatter(
                x=[None], y=[None],
                mode='lines',
                line=dict(color='red', width=2),
                name='Heat Flow',
                showlegend=True
            )
        )
        
        # Add heat flow arrows as individual line traces (no legend)
        for i in range(len(x_sub)):
            # Calculate arrow components
            x_start, y_start = x_sub[i], y_sub[i]
            x_end, y_end = x_start + u_sub[i], y_start + v_sub[i]
            
            # Arrow shaft
            fig.add_trace(
                go.Scatter(
                    x=[x_start, x_end],
                    y=[y_start, y_end],
                    mode='lines',
                    line=dict(color='red', width=2),
                    showlegend=False,
                    hoverinfo='skip'
                )
            )
            
            # Arrow head (simple triangle)
            if abs(u_sub[i]) > 1e-10 or abs(v_sub[i]) > 1e-10:  # Only if vector has magnitude
                # Calculate arrow head
                arrow_length = 0.01
                angle = np.arctan2(v_sub[i], u_sub[i])
                head_angle = 0.5  # radians
                
                # Arrow head points
                x_head1 = x_end - arrow_length * np.cos(angle - head_angle)
                y_head1 = y_end - arrow_length * np.sin(angle - head_angle)
                x_head2 = x_end - arrow_length * np.cos(angle + head_angle)
                y_head2 = y_end - arrow_length * np.sin(angle + head_angle)
                
                fig.add_trace(
                    go.Scatter(
                        x=[x_head1, x_end, x_head2],
                        y=[y_head1, y_end, y_head2],
                        mode='lines',
                        line=dict(color='red', width=2),
                        showlegend=False,
                        hoverinfo='skip'
                    )
                )
        
        fig.update_layout(
            title="Interactive 2D Temperature Distribution with Heat Flow Vectors",
            xaxis_title="X-coordinate",
            yaxis_title="Y-coordinate",
            width=900, height=600,
            legend=dict(
                x=1.02,  # Position legend between plot and colorbar
                y=0.5,
                bgcolor='rgba(255,255,255,0.8)',
                bordercolor='rgba(0,0,0,0.2)',
                borderwidth=1
            )
        )
        
    elif coordinates.shape[1] == 3:  # 3D problem
        x_coords = coordinates[:, 0]
        y_coords = coordinates[:, 1]
        z_coords = coordinates[:, 2]
        
        # Calculate heat flow vectors for 3D
        V = solution.function_space()
        mesh = V.mesh()
        grad_T = project(grad(solution), VectorFunctionSpace(mesh, 'P', 1))
        grad_vals = grad_T.compute_vertex_values(mesh)
        N = coordinates.shape[0]
        grad_x = -grad_vals[:N]    # negative gradient (heat flow x)
        grad_y = -grad_vals[N:2*N] # negative gradient (heat flow y)
        grad_z = -grad_vals[2*N:]  # negative gradient (heat flow z)
        
        # Create 3D interactive plot
        fig = go.Figure()
        
        # Main temperature scatter plot
        fig.add_trace(
            go.Scatter3d(
                x=x_coords,
                y=y_coords,
                z=z_coords,
                mode='markers',
                marker=dict(
                    size=6,
                    color=solution_values,
                    colorscale='Viridis',
                    colorbar=dict(title="Temperature (°C)", x=1.1),
                    opacity=0.7,
                    line=dict(width=1, color='black')
                ),
                name='Temperature',
                hovertemplate='X: %{x:.3f}<br>Y: %{y:.3f}<br>Z: %{z:.3f}<br>Temperature: %{marker.color:.3f}°C<extra></extra>'
            )
        )
        
        # Add heat flow vectors in 3D (heavily subsampled)
        skip = max(1, N//15)  # Even fewer arrows for 3D
        scale_factor = 0.05   # Scale for 3D arrows
        
        # Get subsampled data
        x_sub = x_coords[::skip]
        y_sub = y_coords[::skip]
        z_sub = z_coords[::skip]
        u_sub = grad_x[::skip] * scale_factor
        v_sub = grad_y[::skip] * scale_factor
        w_sub = grad_z[::skip] * scale_factor
        
        # Add a legend entry for heat flow first
        fig.add_trace(
            go.Scatter3d(
                x=[None], y=[None], z=[None],
                mode='lines',
                line=dict(color='red', width=4),
                name='Heat Flow',
                showlegend=True
            )
        )
        
        # Add 3D heat flow arrows (no legend entries)
        for i in range(len(x_sub)):
            if abs(u_sub[i]) > 1e-10 or abs(v_sub[i]) > 1e-10 or abs(w_sub[i]) > 1e-10:
                # Arrow shaft
                fig.add_trace(
                    go.Scatter3d(
                        x=[x_sub[i], x_sub[i] + u_sub[i]],
                        y=[y_sub[i], y_sub[i] + v_sub[i]],
                        z=[z_sub[i], z_sub[i] + w_sub[i]],
                        mode='lines',
                        line=dict(color='red', width=4),
                        showlegend=False,
                        hoverinfo='skip'
                    )
                )
                
                # Simple arrowhead (just a larger point at the end)
                fig.add_trace(
                    go.Scatter3d(
                        x=[x_sub[i] + u_sub[i]],
                        y=[y_sub[i] + v_sub[i]],
                        z=[z_sub[i] + w_sub[i]],
                        mode='markers',
                        marker=dict(color='red', size=3),
                        showlegend=False,
                        hoverinfo='skip'
                    )
                )
        
        fig.update_layout(
            title="Interactive 3D Temperature Distribution with Heat Flow Vectors",
            scene=dict(
                xaxis_title="X-coordinate",
                yaxis_title="Y-coordinate",
                zaxis_title="Z-coordinate",
                camera=dict(
                    eye=dict(x=1.5, y=1.5, z=1.5)  # Better default viewing angle
                )
            ),
            width=900, height=700,
            legend=dict(
                x=0.02,  # Position legend in the top-left corner
                y=0.98,
                bgcolor='rgba(255,255,255,0.8)',
                bordercolor='rgba(0,0,0,0.2)',
                borderwidth=1
            )
        )
    
    # Ensure static directory exists
    if not os.path.exists(STATIC_DIR):
        os.makedirs(STATIC_DIR)

    # Save as interactive HTML file
    plot_filename = os.path.join(STATIC_DIR, 'solution_plot.html')
    fig.write_html(plot_filename)
    
    # Try to save as static PNG for compatibility (optional)
    try:
        png_filename = os.path.join(STATIC_DIR, 'solution_plot.png')
        fig.write_image(png_filename, width=800, height=600)
    except Exception as e:
        print(f"Warning: Could not save static PNG plot: {e}")
        print("Interactive HTML plot is still available.")

    # Also generate the static matplotlib plot
    static_plot_filename = generate_static_solution_plot(solution)
    
    return {
        'interactive': 'solution_plot.html',
        'static': static_plot_filename
    }

def generate_static_solution_plot(solution):
    """Generate a static matplotlib plot (backup function)"""
    
    plt.figure(figsize=(10, 6))
    
    # Get mesh and solution data
    mesh = solution.function_space().mesh()
    coordinates = mesh.coordinates()
    solution_values = solution.compute_vertex_values(mesh)
    
    # Handle different dimensions
    if coordinates.shape[1] == 1:  # 1D problem
        sorted_indices = np.argsort(coordinates[:, 0])
        x_coords = coordinates[sorted_indices, 0]
        y_values = solution_values[sorted_indices]
        plt.plot(x_coords, y_values, 'b-', linewidth=2, marker='o', markersize=4)
        plt.xlabel('X-coordinate')
        plt.ylabel('Temperature (°C)')
        plt.title('Temperature Distribution (1D)')
        plt.grid(True, alpha=0.3)
        # Heat flow vector (negative gradient)
        # Use finite difference for gradient
        grad_T = np.gradient(y_values, x_coords)
        q = -grad_T  # Heat flux direction
        # Plot arrows
        skip = max(1, len(x_coords)//20)
        plt.quiver(x_coords[::skip], y_values[::skip], q[::skip], np.zeros_like(q[::skip]),
                   angles='xy', scale_units='xy', scale=1, color='r', width=0.005)
        plt.legend(['Temperature', 'Heat Flow'])
    elif coordinates.shape[1] == 2:  # 2D problem
        try:
            c = plot(solution, cmap='viridis')
            plt.colorbar(c, label='Temperature (°C)')
            plt.xlabel('X-coordinate')
            plt.ylabel('Y-coordinate')
            plt.title('Temperature Distribution (2D)')
        except:
            x_coords = coordinates[:, 0]
            y_coords = coordinates[:, 1]
            scatter = plt.scatter(x_coords, y_coords, c=solution_values, 
                                cmap='viridis', s=50, edgecolors='black', linewidth=0.5)
            plt.colorbar(scatter, label='Temperature (°C)')
            plt.xlabel('X-coordinate')
            plt.ylabel('Y-coordinate')
            plt.title('Temperature Distribution (2D)')
        # Overlay heat flow vectors
        V = solution.function_space()
        mesh = V.mesh()
        # Interpolate gradient at mesh vertices
        grad_T = project(grad(solution), VectorFunctionSpace(mesh, 'P', 1))
        grad_vals = grad_T.compute_vertex_values(mesh)
        # grad_vals shape: (2*N,) for 2D, where N = number of vertices
        N = coordinates.shape[0]
        grad_x = -grad_vals[:N]  # negative gradient (heat flow x)
        grad_y = -grad_vals[N:]  # negative gradient (heat flow y)
        skip = max(1, N//30)
        plt.quiver(coordinates[::skip, 0], coordinates[::skip, 1], grad_x[::skip], grad_y[::skip],
                   color='r', scale=30, width=0.005)
        plt.legend(['Temperature', 'Heat Flow'])
    elif coordinates.shape[1] == 3:  # 3D problem
        x_coords = coordinates[:, 0]
        y_coords = coordinates[:, 1]
        z_coords = coordinates[:, 2]
        z_middle = (z_coords.min() + z_coords.max()) / 2
        tolerance = (z_coords.max() - z_coords.min()) * 0.1
        mask = np.abs(z_coords - z_middle) < tolerance
        scatter = plt.scatter(x_coords[mask], y_coords[mask], c=solution_values[mask], 
                            cmap='viridis', s=50, edgecolors='black', linewidth=0.5)
        plt.colorbar(scatter, label='Temperature (°C)')
        plt.xlabel('X-coordinate')
        plt.ylabel('Y-coordinate')
        plt.title(f'Temperature Distribution (3D slice at z≈{z_middle:.3f})')
        # Overlay heat flow vectors on slice
        V = solution.function_space()
        mesh = V.mesh()
        grad_T = project(grad(solution), VectorFunctionSpace(mesh, 'P', 1))
        grad_vals = grad_T.compute_vertex_values(mesh)
        N = coordinates.shape[0]
        grad_x = -grad_vals[:N]
        grad_y = -grad_vals[N:2*N]
        grad_z = -grad_vals[2*N:]
        skip = max(1, np.sum(mask)//30)
        plt.quiver(x_coords[mask][::skip], y_coords[mask][::skip], grad_x[mask][::skip], grad_y[mask][::skip],
                   color='r', scale=30, width=0.005)
        plt.legend(['Temperature', 'Heat Flow'])
    
    if not os.path.exists(STATIC_DIR):
        os.makedirs(STATIC_DIR)

    plot_filename = os.path.join(STATIC_DIR, 'solution_plot_static.png')
    plt.savefig(plot_filename, dpi=150, bbox_inches='tight')
    plt.close()

    return 'solution_plot_static.png'