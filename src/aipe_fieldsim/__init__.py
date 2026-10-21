"""aipe_fieldsim: reproducible Gmsh + Elmer FEM field-simulation workflows.

The package wraps three open-source tools (Gmsh for meshing, ElmerGrid for
mesh conversion, ElmerSolver for the finite-element solve), parses their
outputs with pure numpy, and records every case as an auditable directory
(mesh, SIF, logs, numerical results, tool versions, SHA-256 hashes).
"""

__version__ = "0.1.0"

from .tools import ToolInfo, find_gmsh, find_elmer, qualify_all  # noqa: F401
from .elmer import ElmerCase, ElmerTools, ElmerCaseError  # noqa: F401
from .sif import Sif, Section  # noqa: F401
from .vtu import read_vtu, VtuMesh  # noqa: F401
from .savedata import read_save_scalars, read_save_line  # noqa: F401
