"""Shared summary writer for the parametric models."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Optional

from .. import __version__
from ..elmer import ElmerTools, jsonable, now_rfc3339
from ..postprocess import relative_error

MODEL_SUMMARY_SCHEMA = "aipe.field-simulation/model-summary/0.1"


class Model:
    name = "model"
    title = ""
    description = ""

    def __init__(self, workdir: Path | str, params: Optional[dict] = None, tools: Optional[ElmerTools] = None,
                 plots: bool = True):
        self.workdir = Path(workdir)
        self.workdir.mkdir(parents=True, exist_ok=True)
        self.params = dict(self.default_params())
        if params:
            unknown = set(params) - set(self.params)
            if unknown:
                raise KeyError(f"unknown parameters for {self.name}: {sorted(unknown)}")
            self.params.update(params)
        self.tools = tools or ElmerTools.discover()
        self.plots = plots
        self.started_at = now_rfc3339()

    def default_params(self) -> dict:
        raise NotImplementedError

    @staticmethod
    def result(value: float, unit: str, note: str = "") -> dict:
        return {"value": float(value), "unit": unit, "note": note}

    @staticmethod
    def check(name: str, value: float, tolerance: float, description: str = "") -> dict:
        return {"name": name, "value": float(value), "tolerance": tolerance, "passed": bool(abs(value) <= tolerance),
                "description": description}

    @staticmethod
    def rel(value: float, reference: float) -> float:
        return relative_error(value, reference)

    def write_summary(self, status: str, results: dict, checks: list, limitations: list, timing: dict,
                      extra: Optional[dict] = None, evidence: Optional[list] = None, status_note: str = "") -> dict:
        summary = {
            "schema": MODEL_SUMMARY_SCHEMA,
            "package": {"name": "aipe-fieldsim", "version": __version__},
            "model": self.name,
            "title": self.title,
            "description": self.description,
            "status": status,
            "status_note": status_note,
            "created_at": self.started_at,
            "finished_at": now_rfc3339(),
            "parameters": self.params,
            "results": results,
            "checks": checks,
            "all_checks_passed": all(c["passed"] for c in checks) if checks else None,
            "limitations": limitations,
            "timing": timing,
            "tools": self.tools.manifest(),
            "evidence": evidence or [],
        }
        if extra:
            summary.update(extra)
        (self.workdir / "summary.json").write_text(json.dumps(jsonable(summary), indent=2), encoding="utf-8")
        return summary
