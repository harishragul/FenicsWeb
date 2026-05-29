"""Tests for mesh generation — mesh/mesh.py"""
import json
import pytest
from fenics import (
    UnitIntervalMesh, UnitSquareMesh, UnitCubeMesh,
    IntervalMesh, RectangleMesh, BoxMesh, Point,
)


# ---------------------------------------------------------------------------
# Built-in FEniCS mesh constructors
# ---------------------------------------------------------------------------

def test_unit_interval_mesh():
    mesh = UnitIntervalMesh(10)
    assert mesh.num_cells() == 10
    assert mesh.geometry().dim() == 1


def test_unit_square_mesh():
    mesh = UnitSquareMesh(4, 4)
    assert mesh.num_cells() > 0
    assert mesh.geometry().dim() == 2


def test_unit_cube_mesh():
    mesh = UnitCubeMesh(2, 2, 2)
    assert mesh.num_cells() > 0
    assert mesh.geometry().dim() == 3


def test_interval_mesh_custom_domain():
    mesh = IntervalMesh(20, 0.0, 5.0)
    coords = mesh.coordinates()
    assert abs(coords.min() - 0.0) < 1e-12
    assert abs(coords.max() - 5.0) < 1e-12
    assert mesh.num_cells() == 20


def test_rectangle_mesh():
    mesh = RectangleMesh(Point(0, 0), Point(2, 1), 10, 5)
    assert mesh.geometry().dim() == 2
    coords = mesh.coordinates()
    assert coords[:, 0].max() == pytest.approx(2.0)
    assert coords[:, 1].max() == pytest.approx(1.0)


def test_box_mesh():
    mesh = BoxMesh(Point(0, 0, 0), Point(1, 1, 1), 3, 3, 3)
    assert mesh.geometry().dim() == 3
    assert mesh.num_cells() > 0


# ---------------------------------------------------------------------------
# mesh.mesh helper functions
# ---------------------------------------------------------------------------

def test_generate_interactive_mesh_plot_returns_valid_json():
    from mesh.mesh import generate_interactive_mesh_plot
    mesh = UnitSquareMesh(4, 4)
    result = generate_interactive_mesh_plot(mesh)
    data = json.loads(result)
    assert "data" in data


def test_show_mesh_roundtrip_rectangle():
    """serialize → deserialize a rectangle mesh string and check dimensions."""
    from mesh.mesh import show_mesh
    mesh_str = "rectangle,0,0,1,1,4,4"
    mesh = show_mesh(mesh_str)
    assert mesh.geometry().dim() == 2


def test_show_mesh_roundtrip_interval():
    from mesh.mesh import show_mesh
    mesh_str = "interval,10,0.0,1.0"
    mesh = show_mesh(mesh_str)
    assert mesh.geometry().dim() == 1
    assert mesh.num_cells() == 10
