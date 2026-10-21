"""End-to-end tests that execute Gmsh, ElmerGrid and ElmerSolver (marker: solver).

They run the coarsest level of each benchmark plus small instances of the two DAB models and are
skipped automatically when the tools are not found (see conftest.py).
"""
import json

import pytest

pytestmark = pytest.mark.solver


def _comparison(summary, quantity):
    return next(c for c in summary["comparisons"] if c["quantity"] == quantity)


def test_elmer_qualification(solver_tools, tmp_path):
    from aipe_fieldsim.tools import qualify_elmer

    result = qualify_elmer(tmp_path)
    assert result["passed"], result
    record = json.loads((tmp_path / "elmer-qualification" / "case.json").read_text())
    assert record["tools"]["elmersolver"]["version"]
    assert any(a["path"] == "case.sif" for a in record["artifacts"])
    assert record["mesh"]["boundaries"]["4"]["bbox_max"][0] == 0.0  # side 4 is the left edge x = 0


def test_slab_benchmark_coarse(solver_tools, tmp_path):
    from aipe_fieldsim.benchmarks import SlabBenchmark

    summary = SlabBenchmark(tmp_path, tools=solver_tools, plots=False).run(levels=(0,))
    assert _comparison(summary, "T_face_K")["rel_error_finest"] < 1e-8
    assert _comparison(summary, "Q_convective_W_per_m")["rel_error_finest"] < 1e-8
    assert _comparison(summary, "Q_dirichlet_W_per_m")["rel_error_finest"] < 1e-8
    assert summary["all_checks_passed"]
    assert (tmp_path / "summary.json").is_file()


def test_skin_depth_benchmark_coarse(solver_tools, tmp_path):
    from aipe_fieldsim.benchmarks import SkinDepthBenchmark

    summary = SkinDepthBenchmark(tmp_path, tools=solver_tools, plots=False).run(levels=(0,))
    assert _comparison(summary, "loss_per_depth_W_per_m")["rel_error_finest"] < 5e-3
    assert _comparison(summary, "J_surface_A_per_m2")["rel_error_finest"] < 1e-6
    assert summary["all_checks_passed"]


def test_coax_benchmark_coarse(solver_tools, tmp_path):
    from aipe_fieldsim.benchmarks import CoaxBenchmark

    summary = CoaxBenchmark(tmp_path, tools=solver_tools, plots=False).run(levels=(0,))
    assert _comparison(summary, "L_per_length_H_per_m")["rel_error_finest"] < 0.05
    assert summary["all_checks_passed"]


def test_fin_benchmark_coarse(solver_tools, tmp_path):
    from aipe_fieldsim.benchmarks import FinBenchmark

    summary = FinBenchmark(tmp_path, tools=solver_tools, plots=False).run(levels=(0,))
    assert _comparison(summary, "heat_rate_W")["rel_error_finest"] < 0.02
    assert _comparison(summary, "T_tip_K")["rel_error_finest"] < 1e-3
    assert summary["all_checks_passed"]


def test_heatsink_model_small(solver_tools, tmp_path):
    from aipe_fieldsim.models import HeatSinkModel

    model = HeatSinkModel(tmp_path, params={"h_levels_m": [4.0e-3], "fin_count": 4, "base_length_m": 60e-3}, tools=solver_tools, plots=False)
    summary = model.run(levels=(0,))
    assert summary["status"] == "executed"
    checks = {c["name"]: c for c in summary["checks"]}
    assert checks["energy_conservation"]["passed"]
    assert checks["footprint_area"]["passed"]
    assert summary["results"]["Q_convective_W"]["value"] == pytest.approx(40.0, rel=1e-6)
    assert summary["results"]["T_max_K"]["value"] > 300.0


def test_transformer_model_static_quick(solver_tools, tmp_path):
    from aipe_fieldsim.models import TransformerWindowModel

    model = TransformerWindowModel(tmp_path, params={"h0_core_m": 4e-3, "h0_window_m": 2e-3}, tools=solver_tools, plots=False)
    summary = model.run(harmonic=False, levels=(0,))
    assert summary["status"] == "executed"
    checks = {c["name"]: c for c in summary["checks"]}
    assert checks["ampere_turn_balance_right_window"]["passed"]
    assert checks["primary_block_current"]["passed"]
    fem = summary["results"]["L_leakage_primary_H"]["value"]
    est = summary["results"]["L_leakage_primary_1d_estimate_H"]["value"]
    assert abs(fem - est) / est < 0.5  # same order of magnitude as the 1-D MMF estimate
    assert summary["results"]["L_magnetising_primary_H"]["value"] > fem
