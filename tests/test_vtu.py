import math

import numpy as np
import pytest

from aipe_fieldsim.postprocess import (MU0, boundary_integral, cell_measures, convergence_table, integrate_cell,
                                       magnetic_energy_2d, tet_volumes, tri_gradient, triangle_areas)
from aipe_fieldsim.vtu import VtuMesh, read_vtu, write_vtu_appended, write_vtu_ascii


def unit_square_mesh():
    """Unit square: two triangles + one quad on top (0..1 x 0..2), boundary lines with GeometryIds offset 100."""
    points = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 2, 0], [1, 2, 0]], dtype=float)
    # cells: tri(0,1,2), tri(0,2,3), quad(3,2,5,4), line bottom (0,1) bc 1, line top (4,5) bc 2
    connectivity = np.array([0, 1, 2, 0, 2, 3, 3, 2, 5, 4, 0, 1, 4, 5])
    offsets = np.array([3, 6, 10, 12, 14])
    types = np.array([5, 5, 9, 3, 3])
    u = points[:, 0] + 2 * points[:, 1]  # linear field
    cell_data = {"GeometryIds": np.array([1, 1, 2, 101, 102]), "cellval": np.array([1.0, 2.0, 3.0, 0.0, 0.0])}
    return VtuMesh(points, connectivity, offsets, types, {"u": u, "vec": np.column_stack([u, -u, 0 * u])}, cell_data)


@pytest.mark.parametrize("writer", [write_vtu_ascii, write_vtu_appended])
def test_roundtrip(tmp_path, writer):
    mesh = unit_square_mesh()
    path = tmp_path / "m.vtu"
    writer(path, mesh)
    back = read_vtu(path)
    assert back.n_points == 6 and back.n_cells == 5
    np.testing.assert_allclose(back.points, mesh.points)
    np.testing.assert_array_equal(back.connectivity, mesh.connectivity)
    np.testing.assert_array_equal(back.offsets, mesh.offsets)
    np.testing.assert_array_equal(back.types, mesh.types)
    np.testing.assert_allclose(back.point_data["u"], mesh.point_data["u"])
    assert back.point_data["vec"].shape == (6, 3)
    np.testing.assert_array_equal(back.geometry_ids(), mesh.cell_data["GeometryIds"])
    assert back.dimension == 2
    assert back.summary()["cells_by_type"] == {"triangle": 2, "quad": 1, "line": 2}


def test_appended_uint64_header(tmp_path):
    mesh = unit_square_mesh()
    path = tmp_path / "m64.vtu"
    write_vtu_appended(path, mesh, header_type="UInt64")
    back = read_vtu(path)
    np.testing.assert_allclose(back.point_data["u"], mesh.point_data["u"])


def test_measures_and_integrals():
    mesh = unit_square_mesh()
    m = cell_measures(mesh)
    np.testing.assert_allclose(m, [0.5, 0.5, 1.0, 1.0, 1.0])
    assert math.isclose(integrate_cell(mesh, mesh.cell_data["cellval"], mesh.cell_dims() == 2), 0.5 + 1.0 + 3.0)
    integral, length = boundary_integral(mesh, mesh.point_data["u"], 2, offset=100)
    # top edge y=2: u = x + 4, integral over x in [0,1] = 4.5
    assert math.isclose(integral, 4.5) and math.isclose(length, 1.0)
    integral, _ = boundary_integral(mesh, mesh.point_data["u"], 1, offset=100)
    assert math.isclose(integral, 0.5)
    body = mesh.body_of_cells()
    np.testing.assert_array_equal(body, [1, 1, 2, -1, -1])


def test_gradient_and_energy():
    mesh = unit_square_mesh()
    idx, tris = mesh.cells_of_type(5)
    g = tri_gradient(mesh.points, tris, mesh.point_data["u"])
    np.testing.assert_allclose(g, [[1, 2], [1, 2]])
    np.testing.assert_allclose(triangle_areas(mesh.points, tris), [0.5, 0.5])
    # Az = B0 * y gives uniform B_x = B0; energy per depth = B0^2/(2 mu0) * total area (2.0)
    b0 = 0.01
    energy = magnetic_energy_2d(mesh, b0 * mesh.points[:, 1], {1: 1.0, 2: 1.0})
    assert math.isclose(energy["total_J_per_m"], b0 ** 2 / (2 * MU0) * 2.0, rel_tol=1e-12)
    assert math.isclose(energy["per_body_J_per_m"][2], b0 ** 2 / (2 * MU0) * 1.0, rel_tol=1e-12)


def test_tet_volume_and_3d_boundary():
    points = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)
    connectivity = np.array([0, 1, 2, 3, 0, 1, 2])
    mesh = VtuMesh(points, connectivity, np.array([4, 7]), np.array([10, 5]), {"T": np.array([1.0, 1.0, 1.0, 5.0])},
                   {"GeometryIds": np.array([1, 102])})
    assert math.isclose(tet_volumes(points, np.array([[0, 1, 2, 3]]))[0], 1 / 6)
    assert mesh.dimension == 3
    integral, area = boundary_integral(mesh, mesh.point_data["T"], 2, offset=100)
    assert math.isclose(area, 0.5) and math.isclose(integral, 0.5)


def test_convergence_table_orders():
    levels = [{"level": i, "h": 0.1 / 2 ** i, "elements": 10 * 4 ** i, "v": 1.0 + 0.5 * (0.1 / 2 ** i) ** 2} for i in range(3)]
    rows = convergence_table(levels, "v", exact=1.0)
    assert math.isclose(rows[2]["observed_order"], 2.0, rel_tol=1e-9)
    assert rows[0]["rel_error"] > rows[1]["rel_error"] > rows[2]["rel_error"]
    rows = convergence_table(levels, "v")
    assert math.isclose(rows[2]["observed_order"], 2.0, rel_tol=1e-9)
    assert math.isclose(rows[2]["richardson_extrapolated"], 1.0, rel_tol=1e-9)
