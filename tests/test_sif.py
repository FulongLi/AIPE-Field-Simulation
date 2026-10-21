from aipe_fieldsim.sif import (Bare, Matrix, Raw, Section, Sif, calcfields_solver, format_entry, header, heat_solver,
                               mgdyn2d_solver, result_output, save_line, save_scalars, simulation)


def test_format_entry_types():
    assert format_entry("Flag", True) == "Flag = Logical True"
    assert format_entry("Flag", False) == "Flag = Logical False"
    assert format_entry("Count", 3) == "Count = 3"
    assert format_entry("Value", 1.5e-6) == "Value = 1.5e-06"
    assert format_entry("Name", '"copper"') == 'Name = "copper"'
    assert format_entry("Target Boundaries", [1, 2]) == "Target Boundaries(2) = 1 2"
    assert format_entry("Polyline Coordinates", Matrix([[0.0, 0.0], [1.0, 0.5]])) == "Polyline Coordinates(2,2) = 0.0 0.0 1.0 0.5"
    assert format_entry("Density", Raw('Variable Temperature; Real MATC "2*tx"')) == 'Density = Variable Temperature; Real MATC "2*tx"'
    assert format_entry("Mesh DB", Bare('"." "mesh"')) == 'Mesh DB "." "mesh"'


def test_section_render_and_comment():
    s = Section("Boundary Condition", 2, comment="convective face", **{"Target Boundaries": [2], "Heat Transfer Coefficient": 25.0})
    text = s.render()
    assert text.splitlines()[0] == "! convective face"
    assert "Boundary Condition 2" in text
    assert "  Target Boundaries(1) = 2" in text
    assert text.strip().endswith("End")


def test_header_uses_bare_statements():
    text = header("mesh").render()
    assert 'Mesh DB "." "mesh"' in text
    assert "CHECK KEYWORDS Warn" in text
    assert "=" not in text


def test_sif_order_find_count():
    sif = Sif(preamble=["$ I = 10.0"])
    sif.add(header())
    sif.add(simulation("Cartesian 3D", extra={"Angular Frequency": 100.0}))
    sif.add(Section("Body", 1, Equation=1, Material=1))
    sif.add(Section("Body", 2, Equation=1, Material=2))
    text = sif.render()
    assert text.startswith("$ I = 10.0\n\nHeader")
    assert text.index("Simulation") < text.index("Body 1") < text.index("Body 2")
    assert sif.count("Body") == 2
    assert sif.find("Body", 2)["Material"] == 2
    assert "Angular Frequency = 100.0" in text
    assert "Coordinate System = Cartesian 3D" in text


def test_solver_builders():
    heat = heat_solver(1, linear_system="iterative", calculate_loads=True).render()
    assert 'Procedure = "HeatSolve" "HeatSolver"' in heat
    assert "Linear System Iterative Method = BiCGStabl" in heat
    assert "Calculate Loads = Logical True" in heat
    harmonic = mgdyn2d_solver(1, harmonic=True).render()
    assert "MagnetoDynamics2DHarmonic" in harmonic
    assert "Variable = Potential[Potential Re:1 Potential Im:1]" in harmonic
    assert "MagnetoDynamicsCalcFields" in calcfields_solver(2, joule_heating=True).render()
    scalars = save_scalars(3, "scalars.dat", [("Temperature", "max"), ("Temperature", "boundary int", "Conv Side"),
                                             ("Temperature", "diffusive flux", "Conv Side", "Heat Conductivity")]).render()
    assert 'Filename = "scalars.dat"' in scalars
    assert "Mask Name 2 = Conv Side" in scalars
    assert "Coefficient 3 = Heat Conductivity" in scalars
    assert "Mask Name 1" not in scalars
    line = save_line(4, "profile.dat", [(0.0, 0.01), (0.1, 0.01)], 20).render()
    assert "Polyline Coordinates(2,2) = 0.0 0.01 0.1 0.01" in line
    assert "Polyline Divisions(1) = 20" in line
    out = result_output(5, "case", binary=False).render()
    assert "Binary Output = Logical False" in out
    assert "Save Geometry Ids = Logical True" in out
