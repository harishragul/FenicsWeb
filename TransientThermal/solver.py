import os
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.tri as mtri
import plotly.graph_objects as go
import plotly.utils
from fenics import (
    FunctionSpace, TrialFunction, TestFunction, Function,
    Constant, DirichletBC, Expression, interpolate, errornorm,
    dot, grad, dx, solve, lhs, rhs,
    UnitIntervalMesh,
)
from FenicsWeb.settings import BASE_DIR

STATIC_DIR = os.path.join(BASE_DIR, 'static')

# Initial-condition expressions, keyed by dimension and a fixed set of
# user-selectable choices (never built from raw user input — FEniCS
# Expression strings are JIT-compiled as C++).
INITIAL_CONDITIONS = {
    '1D': {
        'sine': "sin(pi*x[0])",
        'gaussian': "exp(-50*pow(x[0] - 0.5, 2))",
        'constant': "1.0",
    },
    '2D': {
        'sine': "sin(pi*x[0]) * sin(pi*x[1])",
        'gaussian': "exp(-50*(pow(x[0] - 0.5, 2) + pow(x[1] - 0.5, 2)))",
        'constant': "1.0",
    },
}


def solve_transient(dimension, mesh, rho, cp, k, f, T, dt, scheme,
                     initial_condition='sine',
                     left_bc=0.0, right_bc=0.0, bottom_bc=0.0, top_bc=0.0):
    """Solve rho*cp*du/dt = k*div(grad(u)) + f from t=0 to t=T.

    `scheme` is 'backward_euler' (1st-order, F uses grad(u)) or
    'crank_nicolson' (2nd-order, F uses grad(0.5*(u + u_n))).

    Dirichlet BCs are applied on every domain boundary using the mesh's own
    coordinate bounds (left/right for 1D; left/right/bottom/top for 2D),
    matching the pattern used by ConvectionDiffusion.solve_convection_diffusion.

    Returns (times, snapshots, coords, u_final):
      - times: [0.0, dt, 2*dt, ..., T]
      - snapshots: list of vertex-value lists (one per time, aligned with coords)
      - coords: mesh.coordinates() as a nested list
      - u_final: the FEniCS Function holding the solution at t=T
    """
    V = FunctionSpace(mesh, "P", 1)
    u = TrialFunction(V)
    v = TestFunction(V)

    ic_expr = INITIAL_CONDITIONS[dimension][initial_condition]
    u_n = interpolate(Expression(ic_expr, degree=4), V)

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

    nsteps = max(1, round(T / dt))
    times = [0.0]
    snapshots = [u_n.compute_vertex_values(mesh).tolist()]

    for step in range(1, nsteps + 1):
        t = step * dt
        if scheme == 'crank_nicolson':
            u_mid = 0.5 * (u + u_n)
            F = (rho * cp / dt * (u - u_n) * v * dx
                 + k * dot(grad(u_mid), grad(v)) * dx
                 - Constant(f) * v * dx)
        else:
            F = (rho * cp / dt * (u - u_n) * v * dx
                 + k * dot(grad(u), grad(v)) * dx
                 - Constant(f) * v * dx)

        a, L = lhs(F), rhs(F)
        u_sol = Function(V)
        solve(a == L, u_sol, bcs)

        times.append(round(t, 10))
        snapshots.append(u_sol.compute_vertex_values(mesh).tolist())
        u_n.assign(u_sol)

    return times, snapshots, coords.tolist(), u_n


def exact_transient_1d(x, t, k=1.0, rho=1.0, cp=1.0):
    """Exact solution u(x,t) = sin(pi x) * exp(-pi^2 k t / (rho cp)) for the
    1D rod with sin(pi x) initial condition and homogeneous Dirichlet BCs."""
    x = np.asarray(x)
    return np.sin(np.pi * x) * np.exp(-np.pi ** 2 * k * t / (rho * cp))


def scheme_comparison_1d(dt_values=(0.02, 0.01, 0.005), T=0.1, k=1.0, rho=1.0, cp=1.0, n=200):
    """Temporal convergence study: solve the 1D rod benchmark (sin(pi x) IC,
    homogeneous Dirichlet, f=0) with both Backward Euler and Crank-Nicolson
    at each dt in dt_values, and compare against the exact solution at t=T.

    Returns a dict with 'dt_values', 'exact_u_center', and one entry per
    scheme containing:
      - 'errors': L2 error at t=T for each dt
      - 'rates':  estimated temporal convergence rate between consecutive
                   dt's (None for the first entry)
      - 'times' / 'u_center': time series of u(0.5, t) at the smallest dt
    """
    mesh = UnitIntervalMesh(n)
    coords = mesh.coordinates().flatten()
    center_idx = int(np.argmin(np.abs(coords - 0.5)))
    smallest_dt = min(dt_values)

    results = {'dt_values': list(dt_values)}
    for scheme in ('backward_euler', 'crank_nicolson'):
        errors = []
        series = None
        for dt in dt_values:
            times, snapshots, _, u_final = solve_transient(
                '1D', mesh, rho=rho, cp=cp, k=k, f=0.0, T=T, dt=dt,
                scheme=scheme, initial_condition='sine',
                left_bc=0.0, right_bc=0.0,
            )
            u_exact_T = Expression(
                "sin(pi*x[0]) * exp(-pi*pi*k*t/(rho*cp))", degree=4,
                k=k, rho=rho, cp=cp, t=T,
            )
            errors.append(errornorm(u_exact_T, u_final, 'L2'))

            if dt == smallest_dt:
                series = {
                    'times': times,
                    'u_center': [snap[center_idx] for snap in snapshots],
                }

        rates = [None]
        for i in range(1, len(errors)):
            rates.append(float(
                np.log(errors[i - 1] / errors[i]) / np.log(dt_values[i - 1] / dt_values[i])
            ))

        results[scheme] = {'errors': errors, 'rates': rates, **series}

    results['exact_u_center'] = exact_transient_1d(
        0.5, np.array(results['backward_euler']['times']), k, rho, cp,
    ).tolist()
    return results


def generate_transient_plotly(times, snapshots, coords):
    """Build a Plotly figure with one frame per timestep and a slider so the
    1D solution u(x,t) can be scrubbed through in the browser without a page
    reload. Returns the figure as a JSON string."""
    coords = np.array(coords).flatten()
    order = np.argsort(coords)
    x_sorted = coords[order].tolist()
    snapshots = np.array(snapshots)
    vmin, vmax = float(snapshots.min()), float(snapshots.max())
    pad = 0.05 * (vmax - vmin) if vmax > vmin else 0.1

    frames = []
    for t, snap in zip(times, snapshots):
        frames.append(go.Frame(
            data=[go.Scatter(x=x_sorted, y=snap[order].tolist(), mode='lines',
                              line=dict(color='firebrick', width=2))],
            name=f"{t:.4g}",
        ))

    fig = go.Figure(
        data=[go.Scatter(x=x_sorted, y=snapshots[0][order].tolist(), mode='lines',
                          line=dict(color='firebrick', width=2))],
        frames=frames,
    )
    fig.update_layout(
        xaxis_title='x', yaxis_title='u(x, t)',
        yaxis=dict(range=[vmin - pad, vmax + pad]),
        sliders=[{
            'active': 0,
            'currentvalue': {'prefix': 'Time t = '},
            'steps': [
                {'args': [[fr.name], {'mode': 'immediate', 'frame': {'duration': 0, 'redraw': True}}],
                 'label': fr.name, 'method': 'animate'}
                for fr in frames
            ],
        }],
        updatemenus=[{
            'type': 'buttons',
            'showactive': False,
            'buttons': [
                {'label': 'Play', 'method': 'animate',
                 'args': [None, {'frame': {'duration': 200, 'redraw': True}, 'fromcurrent': True}]},
                {'label': 'Pause', 'method': 'animate',
                 'args': [[None], {'frame': {'duration': 0, 'redraw': False}, 'mode': 'immediate'}]},
            ],
        }],
    )
    return json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)


def generate_transient_animation(times, snapshots, coords, label='', output_dir=None, fps=5):
    """Generate an animated GIF of the 1D solution u(x,t) evolving over time.
    Returns the output filename."""
    target_dir = output_dir if output_dir is not None else STATIC_DIR
    os.makedirs(target_dir, exist_ok=True)
    suffix = f"_{label}" if label else ""
    filename = os.path.join(target_dir, f'transient_animation{suffix}.gif')

    coords = np.array(coords).flatten()
    order = np.argsort(coords)
    x_sorted = coords[order]
    snapshots = np.array(snapshots)
    vmin, vmax = float(snapshots.min()), float(snapshots.max())
    pad = 0.05 * (vmax - vmin) if vmax > vmin else 0.1

    fig, ax = plt.subplots(figsize=(6, 4))
    line, = ax.plot(x_sorted, snapshots[0][order], 'b-', linewidth=2)
    ax.set_xlabel('x')
    ax.set_ylabel('u(x, t)')
    ax.set_ylim(vmin - pad, vmax + pad)
    ax.grid(True)
    title = ax.set_title(f"t = {times[0]:.4g}")

    def update(frame):
        line.set_ydata(snapshots[frame][order])
        title.set_text(f"t = {times[frame]:.4g}")
        return line, title

    anim = animation.FuncAnimation(fig, update, frames=len(times), blit=False)
    anim.save(filename, writer='pillow', fps=fps)
    plt.close(fig)
    return os.path.basename(filename)


def generate_transient_snapshot_grid(times, snapshots, coords, cells, label='',
                                      output_dir=None, n_panels=4):
    """For 2D problems, plot the temperature field at `n_panels` evenly
    spaced times (including t=0 and t=T) as a row of tripcolor subplots
    sharing a colorbar. Returns the output filename."""
    target_dir = output_dir if output_dir is not None else STATIC_DIR
    os.makedirs(target_dir, exist_ok=True)
    suffix = f"_{label}" if label else ""
    plot_filename = os.path.join(target_dir, f'transient_snapshots{suffix}.png')

    coords = np.array(coords)
    snapshots = np.array(snapshots)
    cells = np.array(cells)
    triang = mtri.Triangulation(coords[:, 0], coords[:, 1], cells)
    vmin, vmax = float(snapshots.min()), float(snapshots.max())

    n_panels = min(n_panels, len(times))
    indices = np.linspace(0, len(times) - 1, n_panels).round().astype(int)

    fig, axes = plt.subplots(1, n_panels, figsize=(4 * n_panels, 4))
    if n_panels == 1:
        axes = [axes]

    tpc = None
    for ax, idx in zip(axes, indices):
        tpc = ax.tripcolor(triang, snapshots[idx], shading='gouraud',
                            vmin=vmin, vmax=vmax, cmap='viridis')
        ax.set_title(f"t = {times[idx]:.4g}")
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_aspect('equal')
    fig.colorbar(tpc, ax=axes, label='u', shrink=0.8)

    title = 'Transient Temperature Field'
    if label:
        title += f' ({label.replace("_", " ").title()})'
    fig.suptitle(title)
    plt.savefig(plot_filename, dpi=100, bbox_inches='tight')
    plt.close(fig)
    return os.path.basename(plot_filename)


def generate_scheme_comparison_plot(results, output_dir=None):
    """Two-panel comparison figure: (left) u(0.5,t) vs t for Backward Euler,
    Crank-Nicolson, and the exact solution at the smallest dt; (right)
    log-log L2 error vs dt for both schemes with O(dt) and O(dt^2) reference
    slopes. Returns the output filename."""
    target_dir = output_dir if output_dir is not None else STATIC_DIR
    os.makedirs(target_dir, exist_ok=True)
    plot_filename = os.path.join(target_dir, 'scheme_comparison.png')

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    be = results['backward_euler']
    cn = results['crank_nicolson']

    ax1.plot(be['times'], results['exact_u_center'], 'k-', label='Exact', linewidth=2)
    ax1.plot(be['times'], be['u_center'], 'r--o', label='Backward Euler', markersize=3)
    ax1.plot(cn['times'], cn['u_center'], 'b-.s', label='Crank-Nicolson', markersize=3)
    ax1.set_xlabel('t')
    ax1.set_ylabel('u(0.5, t)')
    ax1.set_title('Temperature at x=0.5 vs Time')
    ax1.legend()
    ax1.grid(True)

    dt_values = np.array(results['dt_values'], dtype=float)
    ax2.loglog(dt_values, be['errors'], 'r-o', label='Backward Euler')
    ax2.loglog(dt_values, cn['errors'], 'b-s', label='Crank-Nicolson')
    ref1 = be['errors'][0] * (dt_values / dt_values[0])
    ref2 = cn['errors'][0] * (dt_values / dt_values[0]) ** 2
    ax2.loglog(dt_values, ref1, 'k--', alpha=0.5, label='O(dt)')
    ax2.loglog(dt_values, ref2, 'k:', alpha=0.5, label='O(dt²)')
    ax2.set_xlabel('dt')
    ax2.set_ylabel('L2 error at t=T')
    ax2.set_title('Temporal Convergence')
    ax2.legend()
    ax2.grid(True, which='both')

    plt.suptitle('Backward Euler vs Crank-Nicolson — 1D Transient Heat Equation')
    plt.tight_layout()
    plt.savefig(plot_filename, dpi=100, bbox_inches='tight')
    plt.close(fig)
    return os.path.basename(plot_filename)
