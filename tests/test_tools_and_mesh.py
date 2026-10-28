import json
import os
from pathlib import Path

import pytest

from aipe_fieldsim.elmer import jsonable
from aipe_fieldsim.elmermesh import inspect_elmer_mesh, read_names
from aipe_fieldsim.evidence import analytical_evidence, quantity, simulation_evidence
from aipe_fieldsim.grd import BOUNDARY_SIDES, rectangle_grd
from aipe_fieldsim.hashing import artifact_record, hash_tree, sha256_bytes
from aipe_fieldsim.tools import ELMER_BANNER_RE, ELMER_VERSION_RE, ToolInfo, elmer_env, find_executable, tool_manifest


def test_find_executable_env_precedence(tmp_path, monkeypatch):
    exe = tmp_path / "ElmerSolver.exe"
    exe.write_bytes(b"stub")
    monkeypatch.setenv("AIPE_ELMERSOLVER", str(exe))
    path, source = find_executable("ElmerSolver", "AIPE_ELMERSOLVER")
    assert path == exe and source == "env:AIPE_ELMERSOLVER"
    monkeypatch.delenv("AIPE_ELMERSOLVER")
    monkeypatch.setenv("ELMER_HOME", str(tmp_path))
    (tmp_path / "bin").mkdir()
    (tmp_path / "bin" / "ElmerSolver.exe").write_bytes(b"stub")
    (tmp_path / "bin" / "ElmerSolver").write_bytes(b"stub")
    path, source = find_executable("ElmerSolver", "AIPE_ELMERSOLVER")
    assert source == "env:ELMER_HOME" and path.parent == tmp_path / "bin"
    env = elmer_env(path)
    assert env["ELMER_HOME"] == str(tmp_path)
    assert env["PATH"].startswith(str(tmp_path / "bin") + os.pathsep)


def test_version_regex_and_attribution():
    m = ELMER_VERSION_RE.search("MAIN: Version: 9.0 (Rev: Release, Compiled: 2026-01-20)")
    assert m.groups() == ("9.0", "Release", "2026-01-20")
    banner = ELMER_BANNER_RE.search("ELMER SOLVER (v 26.2) STARTED AT: 2026/10/28 14:06:54")
    assert banner.group(1) == "26.2"
    info = ToolInfo("elmersolver", "Elmer", True)
    assert info.attribution() == {"id": "elmersolver", "version": "unversioned", "name": "Elmer"}


def test_tool_manifest_structure():
    manifest = tool_manifest()
    assert set(manifest) >= {"elmersolver", "elmergrid", "gmsh", "python"}
    assert "packages" in manifest["python"]


def test_rectangle_grd_text():
    text = rectangle_grd(0.0, 0.1, 0.0, 0.02, 10, 4)
    assert "Element Divisions 1 = 10" in text and "Element Divisions 2 = 4" in text
    assert "Subcell Limits 1 = 0.0 0.1" in text
    assert BOUNDARY_SIDES[4].startswith("left")
    with pytest.raises(ValueError):
        rectangle_grd(0, 1, 0, 1, 0, 1)


def test_inspect_elmer_mesh(tmp_path):
    (tmp_path / "mesh.header").write_text("4 2 2\n2\n202 2\n303 2\n")
    (tmp_path / "mesh.nodes").write_text("1 -1 0 0 0\n2 -1 1 0 0\n3 -1 1 1 0\n4 -1 0 1 0\n")
    (tmp_path / "mesh.elements").write_text("1 1 303 1 2 3\n2 2 303 1 3 4\n")
    (tmp_path / "mesh.boundary").write_text("1 1 1 0 202 1 2\n2 2 2 0 202 3 4\n")
    (tmp_path / "mesh.names").write_text("! ----- names for bodies -----\n$ lower = 1\n$ upper = 2\n! ----- names for boundaries -----\n$ bottom = 1\n$ top = 2\n")
    info = inspect_elmer_mesh(tmp_path)
    assert info["nodes"] == 4 and info["bulk_elements"] == 2 and info["dimension"] == 2
    assert info["bodies"] == {"1": 1, "2": 1}
    assert info["boundaries"]["2"]["bbox_min"] == [0.0, 1.0, 0.0]
    assert info["names"]["boundaries"]["top"] == 2
    assert read_names(tmp_path / "missing.names") == {"bodies": {}, "boundaries": {}}


def test_hashing_and_jsonable(tmp_path):
    f = tmp_path / "a.txt"
    f.write_bytes(b"hello")
    rec = artifact_record(f, tmp_path)
    assert rec == {"path": "a.txt", "sha256": sha256_bytes(b"hello"), "bytes": 5}
    (tmp_path / "case.json").write_text("{}")
    assert [r["path"] for r in hash_tree(tmp_path)] == ["a.txt"]
    import numpy as np

    out = jsonable({"a": np.float64(1.5), "b": np.array([1, 2]), "c": Path("x"), "d": float("nan"), "e": np.bool_(True)})
    assert out == {"a": 1.5, "b": [1, 2], "c": "x", "d": None, "e": True}
    json.dumps(out)


def test_evidence_fragments():
    q = quantity(40, "W")
    assert q == {"value": 40.0, "unit": "W"}
    ev = simulation_evidence("evidence.x", "desc", "2026-10-21T13:00:00+08:00", [{"id": "elmersolver", "version": "9.0"}],
                             [{"path": "case.json", "sha256": "ab"}], ["assumption.p"])
    assert ev["kind"] == "simulation" and ev["artifacts"][0]["uri"] == "case.json" and ev["review_status"] == "unreviewed"
    an = analytical_evidence("evidence.a", "formula", "2026-10-21T13:00:00+08:00", ["assumption.p"], "textbook")
    assert an["kind"] == "analytical" and an["provenance"]["inputs"] == ["assumption.p"]
    with pytest.raises(ValueError):
        from aipe_fieldsim.evidence import evidence

        evidence("guess", "e", "d", "t")
