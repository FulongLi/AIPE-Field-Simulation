"""Validate aipe.yaml against the AIPE-Registry capability-manifest schema when a sibling checkout exists."""
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
MANIFEST = REPO / "aipe.yaml"


def load_manifest():
    yaml = pytest.importorskip("yaml")
    return yaml.safe_load(MANIFEST.read_text(encoding="utf-8"))


def test_manifest_identity_and_documents_exist():
    m = load_manifest()
    assert m["id"] == "aipe.field-simulation"
    assert m["type"] == "tool" and m["status"] == "experimental" and m["maturity"] == "prototype"
    assert m["license"] == "MIT"
    for ref in m["documentation"] + m["validation"]["evidence"]:
        if not ref.startswith("https://"):
            assert (REPO / ref).exists(), f"referenced file missing: {ref}"
    tool_ids = {t["id"] for t in m["tools"]}
    assert {"gmsh", "elmer", "openfoam"} <= tool_ids
    openfoam = next(t for t in m["tools"] if t["id"] == "openfoam")
    assert openfoam["optional"] is True and openfoam["limitations"]


def test_manifest_against_registry_schema():
    jsonschema = pytest.importorskip("jsonschema")
    candidates = [REPO.parent / "AIPE-Registry" / "schema" / "capability-manifest.schema.json"]
    schema_path = next((c for c in candidates if c.is_file()), None)
    if schema_path is None:
        pytest.skip("AIPE-Registry schema not available next to this repository")
    import json

    schema = json.loads(schema_path.read_text(encoding="utf-8"))
    validator = jsonschema.Draft202012Validator(schema, format_checker=jsonschema.FormatChecker())
    errors = sorted(validator.iter_errors(load_manifest()), key=lambda e: list(e.absolute_path))
    assert not errors, "\n".join(f"{'/'.join(map(str, e.absolute_path))}: {e.message}" for e in errors)
