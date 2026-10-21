"""Closed-form reference solutions used by the benchmarks and the DAB models.

Every function documents its assumptions; the benchmarks compare FEM results against these
with recorded relative errors. All quantities are SI.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable, Sequence

import numpy as np

MU0 = 4e-7 * math.pi


# --------------------------------------------------------------------------------------
# (a) 1-D slab conduction with uniform generation, Dirichlet at x=0 and convection at x=L
# --------------------------------------------------------------------------------------
@dataclass
class SlabSolution:
    k: float
    L: float
    T0: float
    h: float
    T_inf: float
    q_vol: float

    def __post_init__(self):
        k, L, h, q = self.k, self.L, self.h, self.q_vol
        # T(x) = -q x^2/(2k) + C1 x + T0 ; -k T'(L) = h (T(L) - T_inf)
        self.C1 = (q * L + h * q * L ** 2 / (2 * k) - h * (self.T0 - self.T_inf)) / (k + h * L)

    def T(self, x):
        x = np.asarray(x, dtype=float)
        return -self.q_vol * x ** 2 / (2 * self.k) + self.C1 * x + self.T0

    @property
    def T_face(self) -> float:
        return float(self.T(self.L))

    @property
    def x_max(self) -> float:
        if self.q_vol == 0:
            return 0.0 if self.C1 <= 0 else self.L
        return float(np.clip(self.C1 * self.k / self.q_vol, 0.0, self.L))

    @property
    def T_max(self) -> float:
        return float(max(self.T(self.x_max), self.T0, self.T_face))

    @property
    def flux_convective_W_per_m2(self) -> float:
        """Heat flux leaving through the convective face (positive outward)."""
        return self.h * (self.T_face - self.T_inf)

    @property
    def flux_dirichlet_W_per_m2(self) -> float:
        """Heat flux leaving through the Dirichlet face x=0 (positive outward, i.e. in -x direction)."""
        return self.k * self.C1  # -k T'(0) · n with n = -x

    @property
    def generation_W_per_m2(self) -> float:
        return self.q_vol * self.L


# --------------------------------------------------------------------------------------
# (b) coaxial conductor: inner solid conductor radius a with uniform J, return in shell b..c
# --------------------------------------------------------------------------------------
def coax_inductance_per_length(a: float, b: float, c: float) -> dict:
    """Magnetostatic inductance per unit length including both internal inductances.

    Inner conductor (0..a) carries +I uniformly, shell (b..c) carries -I uniformly, field zero for r>c.
    L' = μ0/(2π) [ 1/4 + ln(b/a) + S ], S = [c⁴ ln(c/b) − c²(c²−b²) + (c⁴−b⁴)/4] / (c²−b²)².
    """
    if not (0 < a < b < c):
        raise ValueError("require 0 < a < b < c")
    shell = (c ** 4 * math.log(c / b) - c ** 2 * (c ** 2 - b ** 2) + (c ** 4 - b ** 4) / 4) / (c ** 2 - b ** 2) ** 2
    inner = 0.25
    gap = math.log(b / a)
    total = MU0 / (2 * math.pi) * (inner + gap + shell)
    return {"L_per_length_H_per_m": total, "inner_H_per_m": MU0 / (2 * math.pi) * inner,
            "gap_H_per_m": MU0 / (2 * math.pi) * gap, "shell_H_per_m": MU0 / (2 * math.pi) * shell}


def coax_b_phi(r, a: float, b: float, c: float, current: float) -> np.ndarray:
    r = np.asarray(r, dtype=float)
    out = np.zeros_like(r)
    inner = r <= a
    gap = (r > a) & (r <= b)
    shell = (r > b) & (r <= c)
    out[inner] = MU0 * current * r[inner] / (2 * math.pi * a ** 2)
    out[gap] = MU0 * current / (2 * math.pi * r[gap])
    out[shell] = MU0 * current / (2 * math.pi * r[shell]) * (c ** 2 - r[shell] ** 2) / (c ** 2 - b ** 2)
    return out


# --------------------------------------------------------------------------------------
# (c) AC eddy currents in a slab, planar Az formulation with Dirichlet Az on both faces
# --------------------------------------------------------------------------------------
def skin_depth(frequency: float, conductivity: float, mu_r: float = 1.0) -> float:
    return math.sqrt(2.0 / (2 * math.pi * frequency * MU0 * mu_r * conductivity))


@dataclass
class SlabEddySolution:
    """Slab |x| <= d, Az(-d) = +A0, Az(+d) = -A0 (real), natural BC on y faces, harmonic at frequency f.

    Az(x) = -A0 sinh(k x)/sinh(k d), k = (1+j)/δ ; J_z = -jωσ Az ; H_y = -(1/μ) dAz/dx.
    Time-average loss per unit area in the y-z plane (W/m²):
    P' = (ω²σA0²/2) δ [sinh(2d/δ) − sin(2d/δ)] / [cosh(2d/δ) − cos(2d/δ)].
    """
    A0: float
    d: float
    frequency: float
    conductivity: float
    mu_r: float = 1.0

    def __post_init__(self):
        self.omega = 2 * math.pi * self.frequency
        self.delta = skin_depth(self.frequency, self.conductivity, self.mu_r)
        self.k = (1 + 1j) / self.delta

    def Az(self, x):
        x = np.asarray(x, dtype=float)
        return -self.A0 * np.sinh(self.k * x) / np.sinh(self.k * self.d)

    def Jz(self, x):
        return -1j * self.omega * self.conductivity * self.Az(x)

    def Hy(self, x):
        x = np.asarray(x, dtype=float)
        return self.A0 * self.k * np.cosh(self.k * x) / np.sinh(self.k * self.d) / (MU0 * self.mu_r)

    @property
    def loss_per_area_W_per_m2(self) -> float:
        u = 2 * self.d / self.delta
        return (self.omega ** 2 * self.conductivity * self.A0 ** 2 / 2) * self.delta * \
            (math.sinh(u) - math.sin(u)) / (math.cosh(u) - math.cos(u))

    def loss_per_area_numeric(self, n: int = 200001) -> float:
        x = np.linspace(-self.d, self.d, n)
        j = np.abs(self.Jz(x)) ** 2 / (2 * self.conductivity)
        return float(np.trapezoid(j, x)) if hasattr(np, "trapezoid") else float(np.trapz(j, x))

    def energy_per_area_J_per_m2(self, n: int = 200001) -> float:
        """Time-average magnetic energy per unit y-z area, ∫ ¼ μ |H|² dx."""
        x = np.linspace(-self.d, self.d, n)
        w = 0.25 * MU0 * self.mu_r * np.abs(self.Hy(x)) ** 2
        return float(np.trapezoid(w, x)) if hasattr(np, "trapezoid") else float(np.trapz(w, x))


# --------------------------------------------------------------------------------------
# (d) straight rectangular fin with convective tip (Incropera, Table 3.4 case A)
# --------------------------------------------------------------------------------------
@dataclass
class RectangularFin:
    k: float
    h: float
    thickness: float
    width: float
    length: float
    T_base: float
    T_inf: float

    def __post_init__(self):
        self.P = 2 * (self.width + self.thickness)
        self.A = self.width * self.thickness
        self.m = math.sqrt(self.h * self.P / (self.k * self.A))
        self.theta_b = self.T_base - self.T_inf
        self.hmk = self.h / (self.m * self.k)
        mL = self.m * self.length
        self.denominator = math.cosh(mL) + self.hmk * math.sinh(mL)

    def theta_ratio(self, x):
        x = np.asarray(x, dtype=float)
        arg = self.m * (self.length - x)
        return (np.cosh(arg) + self.hmk * np.sinh(arg)) / self.denominator

    def T(self, x):
        return self.T_inf + self.theta_b * self.theta_ratio(x)

    @property
    def T_tip(self) -> float:
        return float(self.T(self.length))

    @property
    def heat_rate_W(self) -> float:
        mL = self.m * self.length
        return math.sqrt(self.h * self.P * self.k * self.A) * self.theta_b * \
            (math.sinh(mL) + self.hmk * math.cosh(mL)) / self.denominator

    @property
    def efficiency(self) -> float:
        area = self.P * self.length + self.A
        return self.heat_rate_W / (self.h * area * self.theta_b)

    @property
    def biot(self) -> float:
        return self.h * (self.thickness / 2) / self.k


# --------------------------------------------------------------------------------------
# (e) transformer window estimates
# --------------------------------------------------------------------------------------
def leakage_from_mmf_profile(blocks: Sequence[tuple[float, float]], mean_turn_length: float, window_height: float,
                             reference_current: float, n: int = 20001) -> dict:
    """1-D leakage inductance from the MMF diagram across the window width.

    ``blocks`` is the ordered list of ``(width_m, ampere_turns)`` across the window (ampere-turns = 0 for gaps).
    The field is assumed parallel to the layers with H(x) = F(x)/window_height (long-window approximation,
    core μ→∞, no fringing). L_σ = μ0 · MLT / h_w · ∫ (F(x)/I)² dx.
    """
    widths = np.array([b[0] for b in blocks], dtype=float)
    at = np.array([b[1] for b in blocks], dtype=float)
    edges = np.concatenate([[0.0], np.cumsum(widths)])
    x = np.linspace(0.0, edges[-1], n)
    F = np.zeros_like(x)
    cumulative = np.concatenate([[0.0], np.cumsum(at)])
    for i, (w, a) in enumerate(zip(widths, at)):
        sel = (x >= edges[i]) & (x <= edges[i + 1])
        frac = (x[sel] - edges[i]) / w if w > 0 else 0.0
        F[sel] = cumulative[i] + a * frac
    integrand = (F / reference_current) ** 2
    integral = float(np.trapezoid(integrand, x)) if hasattr(np, "trapezoid") else float(np.trapz(integrand, x))
    L = MU0 * mean_turn_length / window_height * integral
    # W' per unit depth = ½ μ0 ∫ H² dA = ½ μ0 ∫ (F/h_w)² h_w dx = ½ μ0 / h_w ∫ F² dx
    energy_per_depth = 0.5 * MU0 / window_height * integral * reference_current ** 2
    return {"L_leakage_H": L, "mmf_closure_A": float(cumulative[-1]),
            "energy_per_unit_depth_J_per_m": energy_per_depth, "window_width_m": float(edges[-1])}


def magnetising_inductance_reluctance(turns: int, mu_r: float, core_area: float, path_length: float,
                                      gap_length: float = 0.0) -> float:
    reluctance = path_length / (MU0 * mu_r * core_area) + gap_length / (MU0 * core_area)
    return turns ** 2 / reluctance


def dowell_factor(xi: float, layers: float) -> float:
    """Dowell AC/DC resistance ratio for a foil layer of normalised thickness xi = h/δ in an m-layer winding."""
    a = (math.sinh(2 * xi) + math.sin(2 * xi)) / (math.cosh(2 * xi) - math.cos(2 * xi))
    b = (math.sinh(xi) - math.sin(xi)) / (math.cosh(xi) + math.cos(xi))
    return xi * (a + 2 * (layers ** 2 - 1) / 3 * b)


def solid_block_ac_factor(width: float, height: float, frequency: float, conductivity: float, mu_r: float = 1.0) -> dict:
    """Skin-depth bookkeeping for a homogenised solid block (reporting aid, not a loss model)."""
    delta = skin_depth(frequency, conductivity, mu_r)
    return {"skin_depth_m": delta, "width_over_delta": width / delta, "height_over_delta": height / delta}
