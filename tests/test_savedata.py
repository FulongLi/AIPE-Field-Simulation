import numpy as np
import pytest

from aipe_fieldsim.savedata import read_names_file, read_save_line, read_save_scalars

NAMES = """Metadata for SaveScalars file: ./scalars.dat
Elmer version: 9.0
Elmer revision: Release
Elmer compilation date: 2026-01-20
Solver input file: case.sif
File started at: 2026/10/21 13:33:00

Variables in columns of matrix:
   1: max: temperature
   2: min: temperature
   3: boundary int: temperature over bc 2
   4: area: temperature over bc 2
   5: boundary sum: temperature loads over bc conv side
   6: sum: temperature loads
   7: res: electromagnetic field energy
"""
DATA = "   3.241529659742E+002   3.000000000000E+002   6.482926829270E+000   2.000000000000E-002  -3.414634146346E+001  -2.000000000002E+003   1.2E-001\n"

LINE_NAMES = """Metadata for SaveLine file: ./heat_line.dat
  bc nodes: 0
  polyline nodes: 3
Data on different columns
  1: Call count
  2: Boundary condition
  3: Node index
  4: coordinate 1
  5: coordinate 2
  6: coordinate 3
  7: temperature
"""
LINE_DATA = """1 1 0   0.0  1.00000000000E-002  0.0  3.00000000000E+002
1 1 1   5.00000000000E-003  1.00000000000E-002  0.0  3.02392502258E+002
1 1 2   1.00000000000E-002  1.00000000000E-002  0.0  3.04661547727E+002
"""


def test_read_names_and_scalars(tmp_path):
    dat = tmp_path / "scalars.dat"
    dat.write_text(DATA)
    (tmp_path / "scalars.dat.names").write_text(NAMES)
    names = read_names_file(tmp_path / "scalars.dat.names")
    assert names[0] == "max: temperature" and names[-1] == "res: electromagnetic field energy"
    table = read_save_scalars(dat)
    assert table.values.shape == (1, 7)
    assert table.last["min: temperature"] == 300.0
    assert table.get("area") == pytest.approx(0.02)
    assert table.match(r"^sum: temperature loads") == pytest.approx(-2000.000000002)
    assert table.match(r"^boundary sum: temperature loads over bc conv side") == pytest.approx(-34.14634146346)
    with pytest.raises(KeyError):
        table.match(r"temperature loads")  # ambiguous
    with pytest.raises(KeyError):
        table.column("temperature")  # ambiguous substring
    assert table.to_dict()["names"] == names


def test_column_count_mismatch(tmp_path):
    dat = tmp_path / "s.dat"
    dat.write_text("1.0 2.0\n")
    (tmp_path / "s.dat.names").write_text(NAMES)
    with pytest.raises(ValueError):
        read_save_scalars(dat)


def test_read_save_line(tmp_path):
    dat = tmp_path / "heat_line.dat"
    dat.write_text(LINE_DATA)
    (tmp_path / "heat_line.dat.names").write_text(LINE_NAMES)
    line = read_save_line(dat)
    assert line.names[-1] == "temperature"
    coords = line.coordinates()
    assert coords.shape == (3, 3)
    np.testing.assert_allclose(coords[:, 0], [0.0, 0.005, 0.01])
    assert line.column("temperature")[0] == 300.0
