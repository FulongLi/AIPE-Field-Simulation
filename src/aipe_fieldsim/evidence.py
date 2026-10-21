"""Helpers producing AIPE-Core Engineering-State style evidence fragments.

Follows ``AIPE-Core/docs/metadata.md``: explicit SI quantities ``{"value": x, "unit": "W"}``,
evidence kinds (``analytical`` / ``simulation`` / ``assumption``), RFC 3339 timestamps of the
actual run, versioned tool attribution ``{id, version, name}``, artifact references with SHA-256,
``review_status: unreviewed``. These fragments are emitted next to the numerical results; a full
Engineering State document is assembled by AIPE-Core tooling, not here.
"""
from __future__ import annotations

from typing import Optional, Sequence


def quantity(value: float, unit: str) -> dict:
    return {"value": float(value), "unit": unit}


def evidence(kind: str, evidence_id: str, description: str, timestamp: str, tools: Sequence[dict] = (),
             artifacts: Sequence[dict] = (), inputs: Sequence[str] = (), source_type: str = "computation",
             source_reference: str = "", locator: Optional[str] = None) -> dict:
    if kind not in {"assumption", "analytical", "simulation", "measurement", "source"}:
        raise ValueError(f"unknown evidence kind {kind}")
    record = {
        "id": evidence_id,
        "kind": kind,
        "description": description,
        "source": {"type": source_type, "reference": source_reference},
        "provenance": {"activity": "aipe_fieldsim run", "tools": list(tools), "inputs": list(inputs)},
        "timestamp": timestamp,
        "artifacts": [{"uri": a["path"], "sha256": a["sha256"]} for a in artifacts],
        "review_status": "unreviewed",
    }
    if locator:
        record["locator"] = locator
    return record


def simulation_evidence(evidence_id: str, description: str, timestamp: str, tools: Sequence[dict],
                        artifacts: Sequence[dict], inputs: Sequence[str] = ()) -> dict:
    return evidence("simulation", evidence_id, description, timestamp, tools, artifacts, inputs,
                    source_type="computation", source_reference="Elmer FEM case directory (case.json)")


def analytical_evidence(evidence_id: str, description: str, timestamp: str, inputs: Sequence[str],
                        reference: str) -> dict:
    return evidence("analytical", evidence_id, description, timestamp, [], [], inputs,
                    source_type="textbook_formula", source_reference=reference)


def assumption_evidence(evidence_id: str, description: str, timestamp: str) -> dict:
    return evidence("assumption", evidence_id, description, timestamp, [], [], [], source_type="assumption",
                    source_reference="declared model assumption")
