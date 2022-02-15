import numpy as np
import pytest

import atlas4py

# -- Fixtures --
@pytest.fixture
def structured_grid():
    grid = atlas4py.StructuredGrid(
        x_spacing=atlas4py.LinearSpacing(-1, 1, 20), y_spacing=atlas4py.LinearSpacing(-1, 1, 21)
    )
    yield grid


@pytest.fixture
def structured_mesh(structured_grid):
    mesh = atlas4py.StructuredMeshGenerator().generate(structured_grid)
    atlas4py.build_edges(mesh)
    atlas4py.build_node_to_edge_connectivity(mesh)
    yield mesh


@pytest.fixture
def structured_function_space(structured_mesh):
    fs = atlas4py.functionspace.CellColumns(structured_mesh)
    yield fs


@pytest.fixture
def structured_in_and_out_fields(structured_function_space):
    fs = structured_function_space
    in_f = fs.create_field(name="my_in_field", levels=1, dtype=np.float64)
    out_f = fs.create_field(name="my_out_field", levels=1, dtype=np.float64)
    yield in_f, out_f


# -- Tests --
def test_version():
    assert isinstance(atlas4py.__version__, str)


def test_grid_generation(structured_grid):
    assert structured_grid.domain.type == "rectangular"
    assert structured_grid.regular == True


def test_mesh_generation(structured_mesh):
    assert structured_mesh.grid.domain.type == "rectangular"


def test_function_space_generation(structured_function_space):
    assert structured_function_space.nb_cells == 380


def test_field_generation(structured_in_and_out_fields):
    in_f, out_f = structured_in_and_out_fields
    out_view = np.asarray(out_f)
    in_view = np.asarray(in_f)
    in_view[...] = np.full_like(in_view, 3.1415)
    out_view[...] = in_view + 1

    assert np.allclose(in_view + 1, out_view)
