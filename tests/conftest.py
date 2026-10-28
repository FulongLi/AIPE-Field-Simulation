"""Shared fixtures. Solver and mesh tests are skipped cleanly when Gmsh or Elmer are not available.

Continuous integration sets ``AIPE_REQUIRE_SOLVER=1`` / ``AIPE_REQUIRE_GMSH=1``: a missing tool then
aborts the run instead of silently skipping the tests that are supposed to execute it.
"""
from __future__ import annotations

import os

import pytest

REQUIRE_ENV = {"solver": "AIPE_REQUIRE_SOLVER", "mesh": "AIPE_REQUIRE_GMSH"}


def _required(marker: str) -> bool:
    return os.environ.get(REQUIRE_ENV[marker], "").strip().lower() in {"1", "true", "yes"}


def _gmsh_reason():
    from aipe_fieldsim.tools import find_gmsh

    info = find_gmsh()
    return None if info.found else "gmsh python module not importable: " + "; ".join(info.notes)


def _tools_or_none():
    try:
        from aipe_fieldsim.elmer import ElmerTools

        tools = ElmerTools.discover()
        reason = _gmsh_reason()
        if reason:
            return None, reason
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
    for marker, unavailable in (("solver", lambda: _tools_or_none()[1]), ("mesh", _gmsh_reason)):
        marked = [item for item in items if marker in item.keywords]
        if not marked:
            continue
        reason = unavailable()
        if reason is None:
            continue
        if _required(marker):
            raise pytest.UsageError(f"{REQUIRE_ENV[marker]} is set but the {marker} tests cannot run: {reason}")
        skip = pytest.mark.skip(reason=f"{marker} tools unavailable: {reason}")
        for item in marked:
            item.add_marker(skip)
