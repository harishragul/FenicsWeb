# Phase 5 — Scientific Data Export & Mesh Import

**Duration:** 2 weeks
**Status:** Not Started
**Depends on:** Phase 1 (steady solver), Phase 2 (convection-diffusion), Phase 3 (transient) — export must cover every solver's output
**Unlocks:** Phase 6 (documentation references these export/import workflows), Phase 7 (reproducibility evidence for JOSS)

---

## Goal

Implement Goals 8 and 9: scientific data export (CSV, VTK, XDMF — all ParaView-compatible) and Gmsh `.msh` mesh import, with reproducible post-processing. These are the two remaining capabilities a researcher needs to use FenicsWeb results outside the browser and to bring real geometries into the solver.

---

## Tasks

### 5.1 CSV Export (Goal 8)

**Purpose:** The simplest, most portable export — degree-of-freedom coordinates paired with solution values, usable in any spreadsheet or plotting tool, and trivial to diff for regression testing.

```python
def export_csv(u_sol, V, filename):
    """Write DOF coordinates and solution values to CSV (x[,y][,z],u)."""
    coords = V.tabulate_dof_coordinates()
    values = u_sol.vector().get_local()
    with open(filename, "w", newline="") as f:
        writer = csv.writer(f)
        dim = coords.shape[1]
        writer.writerow([f"x{i}" for i in range(dim)] + ["u"])
        for c, v in zip(coords, values):
            writer.writerow(list(c) + [v])
```

Add a "Download CSV" link to every results page (1D, 2D, 3D, transient — for transient, one CSV per saved timestep or a single long-format CSV with a `time` column).

---

### 5.2 VTK Export (Goal 8)

VTK (`.pvd` + `.vtu`) is the format ParaView reads natively for static fields.

```python
from fenics import File

def export_vtk(u_sol, filename):
    """Write solution field to legacy VTK (.pvd/.vtu) for ParaView."""
    vtkfile = File(f"{filename}.pvd")
    vtkfile << u_sol
```

Add a "Download VTK" button to every results page.

---

### 5.3 XDMF Export (Goal 8)

XDMF (`.xdmf` + `.h5`) is the preferred modern format — compact, supports time series, and is the recommended FEniCS output for ParaView.

```python
from fenics import XDMFFile

def export_xdmf(u_sol, mesh, filename):
    """Write mesh and solution field to XDMF/HDF5 for ParaView."""
    with XDMFFile(f"{filename}.xdmf") as f:
        f.write(mesh)
        f.write(u_sol)

def export_xdmf_timeseries(filename):
    """Return an open XDMFFile for appending timestepped solutions
    via xdmf.write_checkpoint(u, "u", t, append=True)."""
    return XDMFFile(f"{filename}.xdmf")
```

For the transient solver (Phase 3), write one frame per saved timestep using `write_checkpoint(u, "u", t, append=True)` so ParaView can animate the time series directly.

Add a "Download XDMF" button to every results page.

---

### 5.4 Gmsh Mesh Import (Goal 9)

**Purpose:** Allow solving on real geometries (not just built-in interval/rectangle/box meshes), and demonstrate that the platform can consume standard preprocessing tool output.

The `CustomMesh` form class is already stubbed out (commented) in `mesh/forms.py:62`.

**Implementation steps:**

1. Accept `.msh` (Gmsh format 2 or 4) uploads via `forms.FileField`.
2. Convert with [`meshio`](https://github.com/nschloe/meshio) (pure Python, no external binary required):

   ```python
   import meshio

   def gmsh_to_xdmf(msh_path, xdmf_path):
       """Convert a Gmsh .msh file to XDMF mesh format for FEniCS."""
       mesh = meshio.read(msh_path)
       cell_type = "tetra" if "tetra" in mesh.cells_dict else (
           "triangle" if "triangle" in mesh.cells_dict else "line")
       cells = mesh.cells_dict[cell_type]
       points = mesh.points[:, :2] if cell_type == "triangle" else mesh.points
       out_mesh = meshio.Mesh(points=points, cells={cell_type: cells})
       meshio.write(xdmf_path, out_mesh)
   ```

3. Store the converted `.xdmf`/`.h5` pair in `media/meshes/` and reference it from the `Mesh` model.
4. Load with `XDMFFile(path).read(mesh)` in the solver setup.

**Validation example (required):** Generate a simple Gmsh mesh of the unit square (e.g., via `gmsh` Python API or a checked-in `.msh` file in `examples/meshes/`), import it, and re-run the Phase 4 2D conduction benchmark on the imported mesh. The result must agree with the built-in `UnitSquareMesh` result to within mesh-dependent tolerance — this both validates the import pipeline and demonstrates the workflow end-to-end.

---

## Files to Create

| File | Purpose |
| ---- | ------- |
| `mesh/export.py` | `export_csv`, `export_vtk`, `export_xdmf`, `export_xdmf_timeseries` |
| `mesh/gmsh_import.py` | `gmsh_to_xdmf` conversion via `meshio` |
| `tests/test_export.py` | CSV/VTK/XDMF export correctness tests |
| `tests/test_mesh_import.py` | Gmsh import + re-solve validation test |
| `examples/meshes/unit_square.msh` | Example Gmsh mesh for import validation |
| `docs/usage/export_and_paraview.md` | Export formats + opening results in ParaView |
| `docs/usage/mesh_import.md` | Gmsh import workflow walkthrough |

**Files to modify:**

| File | Change |
| ---- | ------ |
| `mesh/forms.py` | Implement the commented-out `CustomMesh` upload form |
| `mesh/views.py` | Handle `.msh` upload, call `gmsh_to_xdmf`, store result |
| `SteadyStateThermal/views.py`, `ConvectionDiffusion/views.py`, `TransientThermal/views.py` | Wire up export buttons to `mesh/export.py` |
| Solver result templates | Add "Download CSV / VTK / XDMF" links |

---

## Definition of Done

- [ ] CSV export produces `(x[,y][,z],u)` rows for 1D/2D/3D solutions, including transient (long format with `time` column)
- [ ] VTK export opens correctly in ParaView
- [ ] XDMF export opens correctly in ParaView, including a transient time series animation
- [ ] Gmsh `.msh` upload accepted, converted via `meshio`, and solvable
- [ ] Import validation: re-solving the Phase 4 2D benchmark on the imported mesh agrees with the built-in mesh result
- [ ] Export and mesh-import workflows documented in `docs/usage/`
- [ ] `tests/test_export.py` and `tests/test_mesh_import.py` pass in CI
