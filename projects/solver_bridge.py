"""Unpack project config JSON and dispatch to the appropriate solver library."""
import time
from mesh.mesh import show_mesh
from SteadyStateThermal.conduction import solve_conduction


def run_project_solver(project):
    """Run the solver for a project. Return (u_sol, solution_list, solve_time_s)."""
    setup    = project.setup_config
    mesh_cfg = project.mesh_config
    mesh     = show_mesh(mesh_cfg.mesh_str)

    if project.solver_type == 'steady_conduction':
        return _run_steady_conduction(mesh, setup.dimension,
                                      setup.bc_config, setup.solver_params)

    raise ValueError(f"Solver '{project.solver_type}' not yet implemented.")


def _run_steady_conduction(mesh, dimension, bc, sp):
    bc_types       = {face: cfg['type']          for face, cfg in bc.items()}
    neumann_fluxes = {face: cfg.get('value', 0)  for face, cfg in bc.items()
                      if cfg['type'] == 'neumann'}
    robin_hs       = {face: cfg.get('h', 0)      for face, cfg in bc.items()
                      if cfg['type'] == 'robin'}
    robin_u_infs   = {face: cfg.get('u_inf', 0)  for face, cfg in bc.items()
                      if cfg['type'] == 'robin'}

    def _dval(face):
        cfg = bc.get(face, {})
        return float(cfg.get('value', 0.0)) if cfg.get('type') == 'dirichlet' else 0.0

    t0 = time.time()
    u_sol, solution_list = solve_conduction(
        dimension, mesh,
        left_bc=_dval('left'),   right_bc=_dval('right'),
        top_bc=_dval('top'),     bottom_bc=_dval('bottom'),
        front_bc=_dval('front'), back_bc=_dval('back'),
        f=float(sp.get('f', 0.0)),
        k=float(sp.get('k', 1.0)),
        bc_types=bc_types,
        neumann_fluxes=neumann_fluxes,
        robin_hs=robin_hs,
        robin_u_infs=robin_u_infs,
    )
    return u_sol, solution_list, time.time() - t0
