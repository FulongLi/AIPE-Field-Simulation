"""Shared fixtures. Solver tests are skipped cleanly when Gmsh or Elmer are not available."""
from __future__ import annotations

import pytest


def _tools_or_none():
    try:
        from aipe_fieldsim.elmer import ElmerTools
        from aipe_fieldsim.tools import find_gmsh

        tools = ElmerTools.discover()
        if not find_gmsh().found:
            return None, "gmsh python module not importable"
        return tools, None
    except RuntimeError as exc:
        return None, str(exc)


@pytest.fixture(scope="session")
def solver_tools():
    tools, reason = _tools_or_none()
    if tools is None:
        pytest.skip(f"solver tools unavailable: {reason}")
    return tools


def pytest_collection_modifyitems(config, items):
    tools, reason = _tools_or_none()
    if tools is not None:
        return
    skip = pytest.mark.skip(reason=f"solver tools unavailable: {reason}")
    for item in items:
        if "solver" in item.keywords:
            item.add_marker(skip)
