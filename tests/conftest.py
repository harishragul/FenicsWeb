import pytest
from fenics import UnitIntervalMesh, UnitSquareMesh, UnitCubeMesh


@pytest.fixture
def mesh_1d():
    return UnitIntervalMesh(10)


@pytest.fixture
def mesh_2d():
    return UnitSquareMesh(8, 8)


@pytest.fixture
def mesh_3d():
    return UnitCubeMesh(4, 4, 4)
