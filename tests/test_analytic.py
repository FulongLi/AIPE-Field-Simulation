import math

import numpy as np
import pytest

from aipe_fieldsim.analytic import (MU0, RectangularFin, SlabEddySolution, SlabSolution, coax_b_phi,
                                    coax_inductance_per_length, dowell_factor, leakage_from_mmf_profile,
                                    magnetising_inductance_reluctance, skin_depth)


def test_slab_solution_satisfies_bcs_and_energy_balance():
    s = SlabSolution(k=200.0, L=0.1, T0=300.0, h=50.0, T_inf=290.0, q_vol=1.0e6)
    assert s.T(0.0) == pytest.approx(300.0)
    eps = 1e-6
    dTdx_L = (s.T(s.L + eps) - s.T(s.L - eps)) / (2 * eps)  # central difference is exact for a quadratic
    assert -s.k * dTdx_L == pytest.approx(s.h * (s.T_face - s.T_inf), rel=1e-6)
    assert s.generation_W_per_m2 == pytest.approx(s.flux_convective_W_per_m2 + s.flux_dirichlet_W_per_m2, rel=1e-12)
    assert s.T_max >= s.T_face and 0 <= s.x_max <= s.L
    linear = SlabSolution(200.0, 0.1, 300.0, 50.0, 290.0, 0.0)
    assert linear.T_face == pytest.approx(300.0 + (290.0 - 300.0) * 50.0 * 0.1 / (200.0 + 50.0 * 0.1))


def test_coax_inductance_matches_numeric_energy_integral():
    a, b, c, I = 2e-3, 5e-3, 6e-3, 1000.0
    ind = coax_inductance_per_length(a, b, c)
    r = np.linspace(0.0, c, 400001)
    bphi = coax_b_phi(r, a, b, c, I)
    energy = np.trapezoid(bphi ** 2 / (2 * MU0) * 2 * math.pi * r, r)
    assert 2 * energy / I ** 2 == pytest.approx(ind["L_per_length_H_per_m"], rel=1e-6)
    assert ind["inner_H_per_m"] == pytest.approx(MU0 / (8 * math.pi))
    with pytest.raises(ValueError):
        coax_inductance_per_length(5e-3, 2e-3, 6e-3)


def test_skin_depth_and_slab_eddy_closed_form():
    delta = skin_depth(1e5, 5.8e7)
    assert delta == pytest.approx(0.209e-3, rel=2e-3)
    s = SlabEddySolution(A0=1e-6, d=0.5e-3, frequency=1e5, conductivity=5.8e7)
    assert s.loss_per_area_W_per_m2 == pytest.approx(s.loss_per_area_numeric(), rel=1e-7)
    # Ampere: dH_y/dx = J_z
    x = np.array([0.1e-3, 0.3e-3])
    eps = 1e-9
    dh = (s.Hy(x + eps) - s.Hy(x - eps)) / (2 * eps)
    np.testing.assert_allclose(dh, s.Jz(x), rtol=1e-5)
    assert abs(s.Az(np.array([s.d]))[0] + s.A0) < 1e-15
    assert s.energy_per_area_J_per_m2() > 0


def test_rectangular_fin():
    f = RectangularFin(k=200.0, h=25.0, thickness=1e-3, width=40e-3, length=50e-3, T_base=350.0, T_inf=300.0)
    assert f.theta_ratio(0.0) == pytest.approx(1.0)
    eps = 1e-8
    dT = (f.T(eps) - f.T(0.0)) / eps
    assert -f.k * f.A * dT == pytest.approx(f.heat_rate_W, rel=1e-5)
    assert 0 < f.efficiency < 1
    assert f.biot < 1e-3
    assert f.T_tip < f.T_base


def test_leakage_mmf_profile_closed_form():
    # two adjacent blocks, no gaps: L = mu0 N^2 MLT / h_w * (b1 + b2)/3
    N1, I1, N2 = 21, 20.0, 10
    b1 = b2 = 5e-3
    mlt, hw = 0.2, 36e-3
    est = leakage_from_mmf_profile([(b1, N1 * I1), (b2, -N1 * I1)], mlt, hw, I1)
    assert est["L_leakage_H"] == pytest.approx(MU0 * N1 ** 2 * mlt / hw * (b1 + b2) / 3, rel=1e-6)
    assert est["mmf_closure_A"] == pytest.approx(0.0)
    # inserting a gap adds mu0 N^2 MLT g / h_w
    g = 1.5e-3
    est2 = leakage_from_mmf_profile([(b1, N1 * I1), (g, 0.0), (b2, -N1 * I1)], mlt, hw, I1)
    assert est2["L_leakage_H"] - est["L_leakage_H"] == pytest.approx(MU0 * N1 ** 2 * mlt * g / hw, rel=1e-6)
    assert est2["energy_per_unit_depth_J_per_m"] == pytest.approx(0.5 * est2["L_leakage_H"] * I1 ** 2 / mlt, rel=1e-9)


def test_dowell_and_reluctance():
    assert dowell_factor(0.05, 1) == pytest.approx(1.0, abs=1e-4)
    assert dowell_factor(2.0, 3) > dowell_factor(2.0, 1) > 1.0
    L = magnetising_inductance_reluctance(21, 2000.0, 5.4e-4, 0.2)
    assert L == pytest.approx(21 ** 2 * MU0 * 2000.0 * 5.4e-4 / 0.2)
    assert magnetising_inductance_reluctance(21, 2000.0, 5.4e-4, 0.2, gap_length=1e-3) < L
