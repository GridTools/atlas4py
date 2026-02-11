import atlas4py
import numpy as np
import pytest



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
    assert structured_grid.reduced == False
    assert structured_grid.periodic == False
    assert structured_grid.ny == 21
    assert structured_grid.nx == [20 for _ in range(21)]
    assert structured_grid.size == 420
    ll0_0 = structured_grid.lonlat(0, 0)
    ll19_20 = structured_grid.lonlat(19, 20)
    assert [ll0_0.lon, ll0_0.lat] == pytest.approx([-1, -1])
    assert [ll19_20.lon, ll19_20.lat] == pytest.approx([1, 1])
    assert structured_grid.domain.is_global == False


def test_grid_generation_regional():
    grid = atlas4py.StructuredGrid(
        x_spacings=[atlas4py.LinearSpacing(-1, 1, 20) for _ in range(21)],
        y_spacing=atlas4py.LinearSpacing(-1, 1, 21),
    )
    assert grid.regular == True
    assert grid.reduced == False
    assert grid.periodic == False
    assert grid.size == 420
    assert grid.domain.is_global == False


def test_grid_generation_O32():
    grid = atlas4py.Grid("O32")
    with atlas4py.StructuredGrid(grid) as sg:
        assert sg
        assert sg.regular == False
        assert sg.reduced == True
        assert sg.periodic == True
    assert grid.size == 5248
    assert grid.domain.is_global == True


def test_grid_generation_F32():
    grid = atlas4py.Grid("F32")
    with atlas4py.StructuredGrid(grid) as sg:
        assert sg
        assert sg.regular == True
        assert sg.reduced == False
        assert sg.periodic == True
    assert grid.size == 8192
    assert grid.domain.is_global == True


def test_mesh_generation(structured_mesh):
    assert structured_mesh.grid.domain.type == "rectangular"
    assert structured_mesh.cells.size == 380
    assert structured_mesh.nodes.size == 420


def test_mesh_connectivity(structured_mesh):
    mesh = structured_mesh
    assert mesh.cells.node_connectivity.rows == 380
    assert mesh.cells.node_connectivity.cols(0) == 4
    assert mesh.cells.node_connectivity[0, 0] == 20
    assert mesh.cells.node_connectivity[0, 1] == 0
    assert mesh.cells.node_connectivity[0, 2] == 1
    assert mesh.cells.node_connectivity[0, 3] == 21
    assert mesh.cells.node_connectivity[0, 0, 0] == 20
    assert mesh.cells.node_connectivity[0, 0, 1] == 0
    assert mesh.cells.node_connectivity[0, 0, 2] == 1
    assert mesh.cells.node_connectivity[0, 0, 3] == 21

    block = mesh.cells.node_connectivity.block(0)
    assert block.rows == 380
    assert block.cols == 4
    assert block[0, 0] == 20
    assert block[0, 1] == 0
    assert block[0, 2] == 1
    assert block[0, 3] == 21


def test_function_space_generation(structured_function_space):
    assert structured_function_space.nb_cells == 380


def test_field_generation(structured_in_and_out_fields):
    in_f, out_f = structured_in_and_out_fields
    assert in_f.rank == 2
    assert out_f.rank == 2
    assert in_f.shape == [380, 1]
    assert out_f.shape == [380, 1]
    assert np.dtype(in_f.dtype) == np.float64

    out_view = np.asarray(out_f)
    in_view = np.asarray(in_f)
    in_view[...] = np.full_like(in_view, 3.1415)
    out_view[...] = in_view + 1

    assert np.allclose(in_view + 1, out_view)


def test_gmsh_output(structured_mesh):
    import os

    output_file = "test_output.msh"
    with atlas4py.Gmsh(path=output_file) as gmsh:
        gmsh.write(structured_mesh)

    # Check that the file was created and has content
    assert os.path.exists(output_file)

    # Clean up the output file
    os.remove(output_file)

def test_config_from_kwargs():
    config = atlas4py.Config.from_kwargs(option1="value1", option2=42)
    config["option3"] = 3.14
    assert config["option1"] == "value1"
    assert config["option2"] == 42
    assert config["option3"] == 3.14
    assert config.keys == ["option1", "option2", "option3"]

def test_config_from_yaml():
    yaml_string = """
option1: value1
option2: 42
outer:
  inner_option1: 3.14
  inner_option2: true
"""
    config = atlas4py.Config.from_yaml(yaml_string)
    config["outer.inner_option3"] = "added_value"

    assert config.keys == ["option1", "option2", "outer"]
    assert config["option1"] == "value1"
    assert config["option2"] == 42
    assert config["outer.inner_option1"] == 3.14
    assert config["outer.inner_option2"] == True
    outer = config["outer"]
    assert outer["inner_option1"] == 3.14
    assert outer["inner_option2"] == True
    assert outer["inner_option3"] == "added_value"


def test_config_from_file(tmp_path):
    test_config_yaml = tmp_path / "test_config.yaml"
    with open(test_config_yaml, "w") as f:
        f.write(
            """
option1: value1
option2: 42
outer:
    inner_option1: 3.14
    inner_option2: true
"""
        )
    config = atlas4py.Config.from_file(test_config_yaml)
    config["outer.inner_option3"] = "added_value"

    assert config.keys == ["option1", "option2", "outer"]
    assert config["option1"] == "value1"
    assert config["option2"] == 42
    assert config["outer.inner_option1"] == 3.14
    assert config["outer.inner_option2"] == True
    outer = config["outer"]
    assert outer["inner_option1"] == 3.14
    assert outer["inner_option2"] == True
    assert outer["inner_option3"] == "added_value"
