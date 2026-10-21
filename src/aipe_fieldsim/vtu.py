"""A small, dependency-free reader for VTK XML unstructured grids (``.vtu``) as written by Elmer.

Supports ``ascii``, inline ``binary`` (base64, uncompressed) and ``appended`` raw data with
``UInt32``/``UInt64`` block headers, optionally ``vtkZLibDataCompressor`` compressed.
Elmer's ResultOutputSolver writes appended raw, uncompressed, UInt32 headers by default.
"""
from __future__ import annotations

import base64
import re
import xml.etree.ElementTree as ET
import zlib
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

VTK_TYPE_DIM = {1: 0, 3: 1, 21: 1, 5: 2, 9: 2, 22: 2, 23: 2, 8: 2, 10: 3, 12: 3, 13: 3, 14: 3, 24: 3, 25: 3, 11: 3}
VTK_TYPE_NAME = {1: "vertex", 3: "line", 21: "quadratic_line", 5: "triangle", 9: "quad", 22: "quadratic_triangle",
                 23: "quadratic_quad", 10: "tetra", 12: "hexahedron", 13: "wedge", 14: "pyramid",
                 24: "quadratic_tetra", 25: "quadratic_hexahedron"}
_DTYPES = {"Float32": "f4", "Float64": "f8", "Int8": "i1", "UInt8": "u1", "Int16": "i2", "UInt16": "u2",
           "Int32": "i4", "UInt32": "u4", "Int64": "i8", "UInt64": "u8"}


@dataclass
class VtuMesh:
    points: np.ndarray
    connectivity: np.ndarray
    offsets: np.ndarray
    types: np.ndarray
    point_data: dict = field(default_factory=dict)
    cell_data: dict = field(default_factory=dict)
    source: str = ""

    @property
    def n_points(self) -> int:
        return int(self.points.shape[0])

    @property
    def n_cells(self) -> int:
        return int(self.types.shape[0])

    def cell_nodes(self, i: int) -> np.ndarray:
        start = 0 if i == 0 else int(self.offsets[i - 1])
        return self.connectivity[start:int(self.offsets[i])]

    def cell_dims(self) -> np.ndarray:
        lut = np.full(max(VTK_TYPE_DIM) + 1, -1, dtype=int)
        for k, v in VTK_TYPE_DIM.items():
            lut[k] = v
        return lut[self.types]

    @property
    def dimension(self) -> int:
        return int(self.cell_dims().max())

    def cells_of_type(self, vtk_type: int) -> tuple[np.ndarray, np.ndarray]:
        """Return (cell indices, node array of shape (n, nodes_per_cell)) for one VTK type."""
        idx = np.nonzero(self.types == vtk_type)[0]
        if idx.size == 0:
            return idx, np.zeros((0, 0), dtype=int)
        starts = np.where(idx == 0, 0, self.offsets[idx - 1])
        counts = self.offsets[idx] - starts
        n = int(counts[0])
        if not np.all(counts == n):
            raise ValueError("inconsistent node counts for one VTK cell type")
        nodes = np.stack([self.connectivity[s:s + n] for s in starts]) if idx.size else np.zeros((0, n), int)
        return idx, nodes

    def geometry_ids(self) -> np.ndarray:
        if "GeometryIds" not in self.cell_data:
            raise KeyError("GeometryIds not present; enable 'Save Geometry Ids = True' in ResultOutputSolver")
        return self.cell_data["GeometryIds"].astype(int)

    def body_of_cells(self) -> np.ndarray:
        """Body id for bulk cells (cells of the mesh dimension); -1 for boundary cells."""
        dims = self.cell_dims()
        g = self.geometry_ids()
        return np.where(dims == self.dimension, g, -1)

    def boundary_mask(self, bc_id: int, offset: int = 100) -> np.ndarray:
        dims = self.cell_dims()
        g = self.geometry_ids()
        return (dims == self.dimension - 1) & (g == offset + bc_id)

    def summary(self) -> dict:
        dims = self.cell_dims()
        return {
            "points": self.n_points,
            "cells": self.n_cells,
            "cells_by_type": {VTK_TYPE_NAME.get(int(t), str(int(t))): int(n) for t, n in zip(*np.unique(self.types, return_counts=True))},
            "dimension": int(dims.max()),
            "point_data": {k: list(v.shape) for k, v in self.point_data.items()},
            "cell_data": {k: list(v.shape) for k, v in self.cell_data.items()},
        }


def _dtype(type_name: str, byte_order: str) -> np.dtype:
    prefix = "<" if byte_order.lower().startswith("little") else ">"
    return np.dtype(prefix + _DTYPES[type_name])


def _decode_raw_block(raw: bytes, dtype: np.dtype, header_dtype: np.dtype, compressed: bool) -> np.ndarray:
    hs = header_dtype.itemsize
    if not compressed:
        nbytes = int(np.frombuffer(raw[:hs], dtype=header_dtype)[0])
        return np.frombuffer(raw[hs:hs + nbytes], dtype=dtype)
    nblocks = int(np.frombuffer(raw[:hs], dtype=header_dtype)[0])
    header = np.frombuffer(raw[:hs * (3 + nblocks)], dtype=header_dtype)
    sizes = header[3:3 + nblocks].astype(int)
    pos = hs * (3 + nblocks)
    chunks = []
    for size in sizes:
        chunks.append(zlib.decompress(raw[pos:pos + size]))
        pos += size
    return np.frombuffer(b"".join(chunks), dtype=dtype)


def _read_data_array(elem: ET.Element, appended: bytes | None, byte_order: str, header_dtype: np.dtype,
                     compressed: bool) -> np.ndarray:
    dtype = _dtype(elem.get("type", "Float64"), byte_order)
    fmt = elem.get("format", "ascii")
    ncomp = int(elem.get("NumberOfComponents", "1"))
    if fmt == "ascii":
        text = elem.text or ""
        arr = np.array(text.split(), dtype=dtype) if text.strip() else np.zeros(0, dtype=dtype)
    elif fmt == "binary":
        if compressed:
            raise NotImplementedError("inline base64 compressed arrays are not supported")
        raw = base64.b64decode((elem.text or "").strip())
        arr = _decode_raw_block(raw, dtype, header_dtype, False)
    elif fmt == "appended":
        if appended is None:
            raise ValueError("appended data referenced but no AppendedData section found")
        offset = int(elem.get("offset", "0"))
        arr = _decode_raw_block(appended[offset:], dtype, header_dtype, compressed)
    else:
        raise ValueError(f"unknown DataArray format {fmt!r}")
    if ncomp > 1:
        arr = arr.reshape(-1, ncomp)
    return arr


def read_vtu(path: Path | str) -> VtuMesh:
    data = Path(path).read_bytes()
    marker = data.find(b"<AppendedData")
    appended: bytes | None = None
    if marker >= 0:
        xml_bytes = data[:marker] + b"</VTKFile>"
        m = re.search(rb"<AppendedData[^>]*>", data[marker:])
        tag_end = marker + m.end()
        underscore = data.index(b"_", tag_end)
        appended = data[underscore + 1:]
    else:
        xml_bytes = data
    root = ET.fromstring(xml_bytes)
    byte_order = root.get("byte_order", "LittleEndian")
    header_dtype = _dtype(root.get("header_type", "UInt32"), byte_order)
    compressed = root.get("compressor") is not None
    piece = root.find("./UnstructuredGrid/Piece")
    if piece is None:
        raise ValueError("no UnstructuredGrid/Piece in VTU file")

    def read(elem):
        return _read_data_array(elem, appended, byte_order, header_dtype, compressed)

    points = read(piece.find("./Points/DataArray")).astype(float).reshape(-1, 3)
    cells = {e.get("Name"): read(e) for e in piece.findall("./Cells/DataArray")}
    point_data = {e.get("Name"): read(e) for e in piece.findall("./PointData/DataArray")}
    cell_data = {e.get("Name"): read(e) for e in piece.findall("./CellData/DataArray")}
    mesh = VtuMesh(points, cells["connectivity"].astype(int), cells["offsets"].astype(int),
                   cells["types"].astype(int), point_data, cell_data, source=str(path))
    expected = int(piece.get("NumberOfPoints", mesh.n_points))
    if expected != mesh.n_points:
        raise ValueError(f"VTU declares {expected} points but {mesh.n_points} were read")
    return mesh


def write_vtu_ascii(path: Path | str, mesh: VtuMesh) -> None:
    """Minimal ASCII writer (used by the unit tests to round-trip synthetic meshes)."""
    lines = ['<?xml version="1.0"?>', '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">',
             "  <UnstructuredGrid>", f'    <Piece NumberOfPoints="{mesh.n_points}" NumberOfCells="{mesh.n_cells}">']

    def array(name, arr, type_name):
        ncomp = 1 if arr.ndim == 1 else arr.shape[1]
        head = f'<DataArray type="{type_name}" Name="{name}" NumberOfComponents="{ncomp}" format="ascii">' if name else \
            f'<DataArray type="{type_name}" NumberOfComponents="{ncomp}" format="ascii">'
        body = " ".join(repr(float(v)) if type_name.startswith("Float") else str(int(v)) for v in np.asarray(arr).ravel())
        return [f"        {head}", f"          {body}", "        </DataArray>"]

    lines.append("      <PointData>")
    for k, v in mesh.point_data.items():
        lines += array(k, v, "Float64")
    lines.append("      </PointData>")
    lines.append("      <CellData>")
    for k, v in mesh.cell_data.items():
        lines += array(k, v, "Int32" if np.issubdtype(np.asarray(v).dtype, np.integer) else "Float64")
    lines.append("      </CellData>")
    lines.append("      <Points>")
    lines += array(None, mesh.points, "Float64")
    lines.append("      </Points>")
    lines.append("      <Cells>")
    lines += array("connectivity", mesh.connectivity, "Int32")
    lines += array("offsets", mesh.offsets, "Int32")
    lines += array("types", mesh.types, "Int32")
    lines.append("      </Cells>")
    lines += ["    </Piece>", "  </UnstructuredGrid>", "</VTKFile>"]
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_vtu_appended(path: Path | str, mesh: VtuMesh, header_type: str = "UInt32") -> None:
    """Minimal appended-raw writer mirroring Elmer's layout (used by the unit tests)."""
    header_dtype = np.dtype("<" + _DTYPES[header_type])
    blobs: list[bytes] = []
    offsets: list[int] = []

    def register(arr: np.ndarray, dtype: str) -> int:
        payload = np.ascontiguousarray(np.asarray(arr), dtype="<" + _DTYPES[dtype]).tobytes()
        offsets.append(sum(len(b) for b in blobs))
        blobs.append(np.array([len(payload)], dtype=header_dtype).tobytes() + payload)
        return offsets[-1]

    xml = ['<?xml version="1.0"?>',
           f'<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian" header_type="{header_type}">',
           "  <UnstructuredGrid>", f'    <Piece NumberOfPoints="{mesh.n_points}" NumberOfCells="{mesh.n_cells}">',
           "      <PointData>"]
    for k, v in mesh.point_data.items():
        ncomp = 1 if v.ndim == 1 else v.shape[1]
        xml.append(f'        <DataArray type="Float64" Name="{k}" NumberOfComponents="{ncomp}" format="appended" offset="{register(v, "Float64")}"/>')
    xml.append("      </PointData>")
    xml.append("      <CellData>")
    for k, v in mesh.cell_data.items():
        is_int = np.issubdtype(np.asarray(v).dtype, np.integer)
        ncomp = 1 if np.asarray(v).ndim == 1 else v.shape[1]
        t = "Int32" if is_int else "Float64"
        xml.append(f'        <DataArray type="{t}" Name="{k}" NumberOfComponents="{ncomp}" format="appended" offset="{register(v, t)}"/>')
    xml.append("      </CellData>")
    xml.append("      <Points>")
    xml.append(f'        <DataArray type="Float64" NumberOfComponents="3" format="appended" offset="{register(mesh.points, "Float64")}"/>')
    xml.append("      </Points>")
    xml.append("      <Cells>")
    xml.append(f'        <DataArray type="Int32" Name="connectivity" format="appended" offset="{register(mesh.connectivity, "Int32")}"/>')
    xml.append(f'        <DataArray type="Int32" Name="offsets" format="appended" offset="{register(mesh.offsets, "Int32")}"/>')
    xml.append(f'        <DataArray type="Int32" Name="types" format="appended" offset="{register(mesh.types, "Int32")}"/>')
    xml.append("      </Cells>")
    xml += ["    </Piece>", "  </UnstructuredGrid>", '  <AppendedData encoding="raw">', "_"]
    head = ("\n".join(xml)).encode("ascii")
    tail = b"\n  </AppendedData>\n</VTKFile>\n"
    Path(path).write_bytes(head + b"".join(blobs) + tail)
