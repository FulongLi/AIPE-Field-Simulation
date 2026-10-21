"""Common machinery: run refinement levels, compare with analytic values, write summary.json."""
from __future__ import annotations

import json
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Sequence

from .. import __version__
from ..elmer import ElmerTools, jsonable, now_rfc3339
from ..evidence import analytical_evidence, simulation_evidence
from ..hashing import artifact_record
from ..postprocess import convergence_table, relative_error

SUMMARY_SCHEMA = "aipe.field-simulation/benchmark-summary/0.1"


@dataclass
class Level:
    level: int
    h: float
    case_dir: str
    elements: int
    nodes: int
    values: dict
    timing: dict
    warnings: list = field(default_factory=list)
    artifacts: list = field(default_factory=list)


class Benchmark:
    """Subclasses implement ``default_params``, ``analytic``, ``run_level`` and ``comparisons``."""

    name = "benchmark"
    title = ""
    description = ""
    reference = ""

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
        self.levels: list[Level] = []

    # -- to be implemented by subclasses -------------------------------------------------
    def default_params(self) -> dict:
        raise NotImplementedError

    def analytic(self) -> dict:
        """Return {quantity: {"value": v, "unit": u, ...}} computed from closed-form solutions."""
        raise NotImplementedError

    def run_level(self, level: int) -> Level:
        raise NotImplementedError

    def comparisons(self) -> Sequence[tuple[str, str, Optional[str]]]:
        """(fem_key_in_level_values, unit, analytic_key or None)."""
        raise NotImplementedError

    def checks(self, finest: Level) -> list[dict]:
        """Consistency checks on the finest level: [{name, value, tolerance, passed}]."""
        return []

    def finish(self, summary: dict) -> None:
        """Hook for summary plots (profiles, convergence)."""

    # -- driver -------------------------------------------------------------------------
    def run(self, levels: Sequence[int] = (0, 1, 2)) -> dict:
        t0 = time.perf_counter()
        started = now_rfc3339()
        for lvl in levels:
            t_level = time.perf_counter()
            result = self.run_level(int(lvl))
            result.timing.setdefault("level_total_s", round(time.perf_counter() - t_level, 3))
            self.levels.append(result)
        analytic = self.analytic()
        comparisons = []
        for key, unit, akey in self.comparisons():
            exact = analytic[akey]["value"] if akey else None
            rows = convergence_table(
                [{"level": L.level, "h": L.h, "elements": L.elements, "nodes": L.nodes, key: L.values[key]} for L in self.levels],
                key, exact)
            finest = rows[-1]
            entry = {"quantity": key, "unit": unit, "fem_finest": finest["value"], "analytic": exact,
                     "analytic_key": akey, "rel_error_finest": finest.get("rel_error"),
                     "observed_order": finest.get("observed_order"), "levels": rows}
            comparisons.append(entry)
        finest_level = self.levels[-1]
        checks = self.checks(finest_level)
        total = round(time.perf_counter() - t0, 3)
        tool_attribution = self.tools.attributions()
        evidence = [
            analytical_evidence(f"evidence.{self.name}.analytic", f"{self.title}: closed-form reference. {self.reference}",
                                started, [f"assumption.{self.name}.parameters"], self.reference),
        ]
        for L in self.levels:
            evidence.append(simulation_evidence(
                f"evidence.{self.name}.fem.level{L.level}",
                f"{self.title}: Elmer FEM level {L.level}, h = {L.h:.4g} m, {L.elements} bulk elements.",
                started, tool_attribution, L.artifacts, [f"assumption.{self.name}.parameters"]))
        summary = {
            "schema": SUMMARY_SCHEMA,
            "package": {"name": "aipe-fieldsim", "version": __version__},
            "benchmark": self.name,
            "title": self.title,
            "description": self.description,
            "reference": self.reference,
            "created_at": started,
            "finished_at": now_rfc3339(),
            "parameters": self.params,
            "analytic": analytic,
            "levels": [
                {"level": L.level, "h_m": L.h, "elements": L.elements, "nodes": L.nodes, "case_dir": L.case_dir,
                 "values": L.values, "timing": L.timing, "warnings": L.warnings} for L in self.levels
            ],
            "comparisons": comparisons,
            "checks": checks,
            "all_checks_passed": all(c["passed"] for c in checks) if checks else None,
            "timing": {"total_s": total, "per_level_s": [L.timing.get("level_total_s") for L in self.levels]},
            "tools": self.tools.manifest(),
            "evidence": evidence,
        }
        self.finish(summary)
        (self.workdir / "summary.json").write_text(json.dumps(jsonable(summary), indent=2), encoding="utf-8")
        return summary

    # -- helpers for subclasses ---------------------------------------------------------
    @staticmethod
    def artifacts_of(case) -> list[dict]:
        out = []
        for name in ("case.json", "case.sif", "results.json", "elmersolver.log"):
            p = case.case_dir / name
            if p.is_file():
                out.append(artifact_record(p, case.case_dir.parent))
        return out

    @staticmethod
    def check(name: str, value: float, tolerance: float, description: str = "") -> dict:
        return {"name": name, "value": float(value), "tolerance": tolerance, "passed": bool(abs(value) <= tolerance),
                "description": description}

    @staticmethod
    def rel(value: float, reference: float) -> float:
        return relative_error(value, reference)


def scalar(table, operator: str, variable: str, mask: Optional[str] = None, bc_id: Optional[int] = None) -> float:
    """Look up a SaveScalars column robustly.

    Elmer names masked columns either ``<op>: <var> over bc <mask name>`` or ``<op>: <var> over bc <id>``
    (and ``<op>: <var> mask <mask name>`` for body masks), depending on the operator; accept all variants.
    """
    import re as _re

    var = _re.escape(variable.lower())
    op = _re.escape(operator.lower())
    if mask is None:
        return table.match(rf"^{op}: {var}$")
    alternatives = [_re.escape(mask.lower())]
    if bc_id is not None:
        alternatives.append(str(int(bc_id)))
    alt = "|".join(alternatives)
    return table.match(rf"^{op}: {var} (over bc|mask) ({alt})$")
