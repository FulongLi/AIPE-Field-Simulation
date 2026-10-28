"""Solver-independent layers: materials, physics description, its Elmer rendering and the heat-sink geometry."""
import pytest

from aipe_fieldsim.geometry import finned_heatsink_areas
from aipe_fieldsim.materials import MATERIALS, get_material
from aipe_fieldsim.models.heatsink import HEATSINK_DEFAULTS
from aipe_fieldsim.physics import Adiabatic, Convection, HeatFlux, SteadyHeatConduction
from aipe_fieldsim.sif import Sif, heat_conduction_sections, material_section

IDS = {"footprint_1": 1, "footprint_2": 2, "convective": 3, "bottom": 4}


def _problem(**overrides):
    al = get_material("aluminium").with_(**overrides)
    return SteadyHeatConduction({"heatsink": al}, [
        HeatFlux("footprint_1", ("footprint_1",), 1.0e4, note="20 W"),
        HeatFlux("footprint_2", ("footprint_2",), 1.0e4),
        Convection("convective", ("convective",), 25.0, 300.0),
        Adiabatic("bottom_adiabatic", ("bottom",)),
    ])


def test_material_library():
    al = get_material("aluminium")
    assert al.thermal_conductivity_W_per_m_K == 200.0 and al.density_kg_per_m3 == 2700.0
    assert al.with_(thermal_conductivity_W_per_m_K=167.0).thermal_conductivity_W_per_m_K == 167.0
    assert al.thermal_conductivity_W_per_m_K == 200.0  # frozen: with_ returns a copy
    assert set(MATERIALS) >= {"aluminium", "copper"}
    with pytest.raises(KeyError):
        get_material("unobtainium")
    em = material_section(2, get_material("copper"), physics="electromagnetic").render()
    assert "Electric Conductivity = 58000000.0" in em and "Heat Conductivity" not in em


def test_problem_check_and_description():
    problem = _problem()
    problem.check(["heatsink"], IDS)
    desc = problem.describe()
    assert desc["boundary_conditions"][2] == {"type": "Convection", "name": "convective", "boundaries": ["convective"],
                                              "h_W_per_m2_K": 25.0, "T_ambient_K": 300.0, "note": ""}
    with pytest.raises(ValueError, match="unknown boundary"):
        problem.check(["heatsink"], {"footprint_1": 1})
    with pytest.raises(ValueError, match="bodies not in the mesh"):
        problem.check(["base"], IDS)
    twice = SteadyHeatConduction(problem.bodies, problem.boundary_conditions + [Adiabatic("again", ("bottom",))])
    with pytest.raises(ValueError, match="two conditions"):
        twice.check(["heatsink"], IDS)


def test_heat_conduction_sections_render_explicit_sif():
    sif = Sif()
    for s in heat_conduction_sections(_problem(), {"heatsink": 1}, IDS, masks={"convective": "Conv Side"}):
        sif.add(s)
    text = sif.render()
    assert sif.count("Boundary Condition") == 4 and sif.count("Material") == 1
    assert "Target Bodies(1) = 1" in text and "Heat Conductivity = 200.0" in text
    bc3 = sif.find("Boundary Condition", 3).render()
    assert "Heat Transfer Coefficient = 25.0" in bc3 and "External Temperature = 300.0" in bc3
    assert "Conv Side = Logical True" in bc3
    assert sif.find("Boundary Condition", 4).entries == {"Name": '"bottom_adiabatic"', "Target Boundaries": [4]}
    assert "! 20 W" in sif.find("Boundary Condition", 1).render()


def test_heatsink_closed_form_areas_match_recorded_run():
    areas = finned_heatsink_areas(HEATSINK_DEFAULTS)
    assert areas["convective_m2"] == pytest.approx(0.03272, rel=1e-12)  # results/heatsink_h25/summary.json
    assert areas["footprint_m2"] == pytest.approx(16e-3 * 21e-3)


@pytest.mark.mesh
def test_heatsink_geometry_groups_areas_and_mesh(tmp_path):
    from aipe_fieldsim.geometry import finned_heatsink
    from aipe_fieldsim.meshing import generate, gmsh_session, set_mesh_size, write_msh2

    p = dict(HEATSINK_DEFAULTS)
    expected = finned_heatsink_areas(p)
    with gmsh_session("heatsink-test") as gmsh:
        groups = finned_heatsink(gmsh, p)
        assert groups == {"bodies": {"heatsink": 1}, "boundaries": IDS}

        def mass(dim, tag):
            return sum(gmsh.model.occ.getMass(dim, e) for e in gmsh.model.getEntitiesForPhysicalGroup(dim, tag))

        assert mass(3, 1) == pytest.approx(expected["volume_m3"], rel=1e-9)
        for name in ("footprint_1", "footprint_2"):
            assert mass(2, IDS[name]) == pytest.approx(expected["footprint_m2"], rel=1e-9)
        assert mass(2, IDS["convective"]) == pytest.approx(expected["convective_m2"], rel=1e-9)
        assert mass(2, IDS["bottom"]) == pytest.approx(expected["bottom_uncovered_m2"], rel=1e-9)
        set_mesh_size(gmsh, 6e-3, 3e-3)
        stats = generate(gmsh, 3, optimize=False)
        path = write_msh2(gmsh, tmp_path / "heatsink.msh")
    assert stats["elements_by_dim"]["3"] > 1000
    text = path.read_text()
    assert text.startswith("$MeshFormat\n2.2 0 8") and '2 1 "footprint_1"' in text and '3 1 "heatsink"' in text
