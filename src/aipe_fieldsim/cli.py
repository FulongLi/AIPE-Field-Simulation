"""Command line interface: ``aipe-fieldsim``."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from . import __version__
from .elmer import jsonable


def _parse_params(items):
    params = {}
    for item in items or []:
        if "=" not in item:
            raise SystemExit(f"--param expects key=value, got {item!r}")
        key, value = item.split("=", 1)
        try:
            params[key] = json.loads(value)
        except json.JSONDecodeError:
            params[key] = value
    return params


def cmd_doctor(args) -> int:
    from .tools import qualify_all, tool_manifest

    if args.qualify:
        manifest = qualify_all(Path(args.qualify))
    else:
        manifest = tool_manifest()
    if args.json:
        print(json.dumps(jsonable(manifest), indent=2))
    else:
        for key in ("elmersolver", "elmergrid", "gmsh"):
            info = manifest[key]
            state = "found" if info["found"] else "MISSING"
            print(f"{info['name']:<28} {state:<8} version={info.get('version')} source={info.get('source')} path={info.get('path')}")
            for note in info.get("notes", []):
                print(f"    note: {note}")
        py = manifest["python"]
        print(f"{'Python':<28} found    version={py['version']} packages={py['packages']}")
        if "qualification" in manifest:
            for key, res in manifest["qualification"].items():
                print(f"qualification {key}: {'PASSED' if res.get('passed') else 'FAILED'} - {res.get('test')} {res.get('error', '')}")
        print("OpenFOAM (CFD): not integrated - see docs/jobs/openfoam-heatsink-cfd-job.md for the blocked job specification")
    ok = manifest["elmersolver"]["found"] and manifest["elmergrid"]["found"] and manifest["gmsh"]["found"]
    if "qualification" in manifest:
        ok = ok and all(r.get("passed") for r in manifest["qualification"].values())
    return 0 if ok else 1


def cmd_benchmark(args) -> int:
    from .benchmarks import ALIASES, REGISTRY, get_benchmark
    from .elmer import ElmerTools

    names = list(REGISTRY) if args.name == "all" else [ALIASES.get(args.name, args.name)]
    tools = ElmerTools.discover()
    params = _parse_params(args.param)
    levels = list(range(args.levels))
    rc = 0
    for name in names:
        cls = get_benchmark(name)
        bench = cls(Path(args.out) / name, params=params if len(names) == 1 else None, tools=tools, plots=not args.no_plots)
        print(f"[{name}] running levels {levels} in {bench.workdir} ...", flush=True)
        summary = bench.run(levels)
        for comp in summary["comparisons"]:
            err = comp["rel_error_finest"]
            err_text = f"{100 * err:.3g} %" if err is not None else "(no analytic value)"
            print(f"  {comp['quantity']:<40} analytic={comp['analytic']!s:<24} fem={comp['fem_finest']:<24.10g} error={err_text}")
        for ch in summary["checks"]:
            print(f"  check {ch['name']:<40} {ch['value']:.3e} (tol {ch['tolerance']:.1e}) {'ok' if ch['passed'] else 'FAILED'}")
            rc = rc or (0 if ch["passed"] else 2)
        print(f"  total {summary['timing']['total_s']} s; summary: {bench.workdir / 'summary.json'}")
    return rc


def cmd_transformer(args) -> int:
    from .models.transformer_window import TransformerWindowModel

    params = _parse_params(args.param)
    if args.interleave:
        params["interleaved"] = True
    model = TransformerWindowModel(Path(args.out), params=params, plots=not args.no_plots)
    summary = model.run(harmonic=not args.no_harmonic, levels=list(range(args.levels)))
    print(json.dumps(jsonable({k: summary[k] for k in ("status", "results", "checks", "timing") if k in summary}), indent=2))
    return 0 if summary.get("status") == "executed" else 1


def cmd_heatsink(args) -> int:
    from .models.heatsink import HeatSinkModel

    params = _parse_params(args.param)
    params["h_W_per_m2_K"] = args.h
    model = HeatSinkModel(Path(args.out), params=params, plots=not args.no_plots)
    summary = model.run(levels=list(range(args.levels)))
    print(json.dumps(jsonable({k: summary[k] for k in ("status", "results", "checks", "timing") if k in summary}), indent=2))
    return 0 if summary.get("status") == "executed" and summary.get("all_checks_passed", True) else 1


def cmd_report(args) -> int:
    from .report import write_results

    path = write_results(Path(args.runs), Path(args.out))
    print(f"wrote {path}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="aipe-fieldsim", description="Gmsh + Elmer FEM field-simulation workflows with analytic benchmarks")
    p.add_argument("--version", action="version", version=f"aipe-fieldsim {__version__}")
    sub = p.add_subparsers(dest="command", required=True)

    d = sub.add_parser("doctor", help="locate Gmsh/ElmerSolver/ElmerGrid and report versions; --qualify runs tiny known-answer cases")
    d.add_argument("--qualify", metavar="DIR", help="directory for the qualification cases")
    d.add_argument("--json", action="store_true")
    d.set_defaults(func=cmd_doctor)

    b = sub.add_parser("benchmark", help="run an analytic benchmark with a mesh-convergence study")
    b.add_argument("name", help="slab | coax | skin | fin | all (or full names)")
    b.add_argument("--out", default="runs/benchmarks")
    b.add_argument("--levels", type=int, default=3, help="number of refinement levels (default 3)")
    b.add_argument("--param", action="append", help="override a parameter, key=value (JSON value), single benchmark only")
    b.add_argument("--no-plots", action="store_true")
    b.set_defaults(func=cmd_benchmark)

    t = sub.add_parser("transformer", help="2-D E-core transformer window: leakage/magnetising inductance and homogenised AC loss factor")
    t.add_argument("--out", default="runs/transformer")
    t.add_argument("--interleave", action="store_true", help="P-S-P interleaved winding arrangement")
    t.add_argument("--no-harmonic", action="store_true", help="skip the 100 kHz harmonic solve")
    t.add_argument("--levels", type=int, default=3, help="mesh refinement levels for the magnetostatic solves (default 3)")
    t.add_argument("--param", action="append", help="override a parameter, key=value")
    t.add_argument("--no-plots", action="store_true")
    t.set_defaults(func=cmd_transformer)

    hs = sub.add_parser("heatsink", help="3-D finned heat sink with two TO-247 sources and a convective film coefficient")
    hs.add_argument("--out", default="runs/heatsink")
    hs.add_argument("--h", type=float, default=25.0, help="film coefficient W/m^2/K (25 natural, 80 forced)")
    hs.add_argument("--levels", type=int, default=3, help="mesh refinement levels (default 3: 3 mm, 2 mm, 1.4 mm)")
    hs.add_argument("--param", action="append", help="override a parameter, key=value")
    hs.add_argument("--no-plots", action="store_true")
    hs.set_defaults(func=cmd_heatsink)

    r = sub.add_parser("report", help="collect summary.json files into results/benchmarks.md")
    r.add_argument("--runs", default="runs")
    r.add_argument("--out", default="results")
    r.set_defaults(func=cmd_report)
    return p


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
