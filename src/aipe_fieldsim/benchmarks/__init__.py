"""Analytically verifiable benchmarks executed with Gmsh/ElmerGrid/ElmerSolver."""
from __future__ import annotations

from .base import Benchmark  # noqa: F401
from .slab import SlabBenchmark
from .coax import CoaxBenchmark
from .skin import SkinDepthBenchmark
from .fin import FinBenchmark

REGISTRY = {
    SlabBenchmark.name: SlabBenchmark,
    CoaxBenchmark.name: CoaxBenchmark,
    SkinDepthBenchmark.name: SkinDepthBenchmark,
    FinBenchmark.name: FinBenchmark,
}
ALIASES = {"slab": SlabBenchmark.name, "coax": CoaxBenchmark.name, "skin": SkinDepthBenchmark.name, "fin": FinBenchmark.name}


def get_benchmark(name: str):
    key = ALIASES.get(name, name)
    if key not in REGISTRY:
        raise KeyError(f"unknown benchmark {name!r}; choose from {sorted(REGISTRY) + sorted(ALIASES)}")
    return REGISTRY[key]
