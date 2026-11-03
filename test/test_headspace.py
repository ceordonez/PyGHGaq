"""Tests for the headspace functions.

Reference values come from pyheadspace (commit 3bf8d73) run on its
template/DataTemplate.xlsx with sal_psu = 0, so that its salinity
correction does not apply. Volumes: headspace 35 mL, total 140 mL.

pyheadspace uses R = 8.314 J mol-1 K-1, pyghgaq the exact 8.314462618: with
R = 8.314 the results must match pyheadspace exactly (same method); with the
exact R they differ by less than 1e-4 (relative).
"""

import numpy as np
import pytest

import pyghgaq as pyg
import pyghgaq.main
from pyghgaq.functions.read import read_constant
from pyghgaq.functions.units import Q_

PYHEADSPACE_R = 8.314
R_TOLERANCE = 1e-4  # relative difference from using the exact R instead of 8.314

TEMP_EQ = np.array([27.6, 27.1, 27.1, 26.1, 26.1, 27.5])
TEMP_AIR = np.array([28.0, 29.6, 29.6, 26.2, 26.2, 27.0])
PRESS = np.array([1011.0, 1015.2, 1015.2, 1014.7, 1014.7, 1014.3])
VOL_HS = 35
VOL_TOTAL = 140

CH4_HS = np.array([3.68486330579851, 5.87128365867549, 6.24102518371642,
                   6.85200036157363, 7.05921633296716, 6.29175452428756])
CH4_AIR = np.array([2.0, 1.84811408894594, 1.84811408894594,
                    1.84811408894594, 1.87805932886376, 1.84843244508352])
D13CH4_HS = np.array([-63.3328901428707, -65.3765385132401, -65.2598705399452,
                      -66.7065656241243, -66.6431640379154, -66.3041943577163])
D13CH4_AIR = np.array([-50.8697534026854, -50.9366112754057, -50.9366112754057,
                       -50.9366112754057, -50.9957948441973, -50.6348023478947])

CO2_HS = np.array([602.032512820747, 1467.1866426275, 1509.43417602239,
                   1506.70011505771, 1505.95272591095, 1590.46996931891])
CO2_AIR = np.array([430.0, 403.446759746519, 403.446759746519,
                    403.446759746519, 409.154545359643, 404.601416205804])
D13CO2_HS = np.array([-10.7759519982567, -13.93105542263, -13.5225266789104,
                      -14.0141342332687, -13.6436090121443, -13.9550314154478])
D13CO2_AIR = np.array([-8.6548083747887, -8.65002015040175, -8.65002015040175,
                       -8.65002015040175, -8.8195784934138, -8.6885596798752])

REF = {
    "CH4": {
        "c": [0.02770131288012934, 0.06275479364418075, 0.06827144592748455,
              0.07756118444500827, 0.08025936103358021, 0.06856818159211989],
        "d13c": [-75.44098869796245, -71.0924594017186, -70.47150379399436,
                 -71.8121046011868, -71.61804818071116, -72.02592110609297],
        "inputs": (CH4_HS, CH4_AIR, D13CH4_HS, D13CH4_AIR),
    },
    "CO2": {
        "c": [21.06555855794637, 60.935407479408745, 62.84622280740862,
              63.985665711399726, 63.8736360701369, 65.83034680576739],
        "d13c": [-11.358711791346066, -14.401105516300806, -13.943030787885835,
                 -14.47375611850854, -14.063539809541336, -14.393552438232593],
        "inputs": (CO2_HS, CO2_AIR, D13CO2_HS, D13CO2_AIR),
    },
}


@pytest.fixture(params=["pyheadspace R", "exact R"])
def gas_constant(request, monkeypatch):
    """Tolerance for the comparison with pyheadspace, for each value of R."""
    if request.param == "pyheadspace R":
        constants = {**read_constant(), "R": PYHEADSPACE_R}
        monkeypatch.setattr(pyghgaq.main, "read_constant", lambda *a, **k: constants)
        return 1e-10
    return R_TOLERANCE


@pytest.mark.parametrize("gas", ["CH4", "CO2"])
def test_headspace_concentration_matches_pyheadspace(gas, gas_constant):
    c_hs, c_air, _, _ = REF[gas]["inputs"]
    result = pyg.headspace_concentration(
        gas, c_hs, c_air, TEMP_EQ, TEMP_AIR, PRESS, VOL_HS, VOL_TOTAL
    )
    assert result.units == Q_(1, "mmol m^-3").units
    np.testing.assert_allclose(result.magnitude, REF[gas]["c"], rtol=gas_constant)


@pytest.mark.parametrize("gas", ["CH4", "CO2"])
def test_headspace_d13c_matches_pyheadspace(gas, gas_constant):
    c_hs, c_air, d_hs, d_air = REF[gas]["inputs"]
    result = pyg.headspace_d13c(
        gas, c_hs, c_air, d_hs, d_air, TEMP_EQ, TEMP_AIR, PRESS, VOL_HS, VOL_TOTAL
    )
    assert result.units == Q_(1, "permille").units
    np.testing.assert_allclose(result.magnitude, REF[gas]["d13c"], rtol=gas_constant)


def test_gas_constant_is_exact():
    assert read_constant()["R"] == 8.314462618


def test_headspace_quantity_inputs_are_converted():
    plain = pyg.headspace_concentration(
        "CH4", CH4_HS[0], CH4_AIR[0], TEMP_EQ[0], TEMP_AIR[0], PRESS[0],
        VOL_HS, VOL_TOTAL,
    )
    quantities = pyg.headspace_concentration(
        "CH4",
        Q_(CH4_HS[0] * 1000, "ppb"),
        Q_(CH4_AIR[0] * 1000, "ppb"),
        Q_(TEMP_EQ[0] + 273.15, "K"),
        Q_(TEMP_AIR[0] + 273.15, "K"),
        Q_(PRESS[0] * 100, "Pa"),
        Q_(0.035, "L"),
        Q_(0.140, "L"),
    )
    np.testing.assert_allclose(quantities.magnitude, plain.magnitude, rtol=1e-10)


def test_headspace_mass_balance_zero_when_no_gas_added():
    # Gas-free water: after shaking, the gas initially in the headspace is
    # split between headspace and water, so c_hs < c_air and the water
    # concentration before shaking must come out as zero.
    c_air = 2.0
    hcp = pyg.henry_coefficient("CH4", 20).to("mol m^-3 Pa^-1").magnitude
    v_hs, v_w = 35e-6, 105e-6
    rt = 8.314462618 * 293.15
    # n_air = n_hs + n_aq  ->  x_air V_hs / RT = x_hs (V_hs / RT + H V_w)
    c_hs = c_air * (v_hs / rt) / (v_hs / rt + hcp * v_w)
    result = pyg.headspace_concentration("CH4", c_hs, c_air, 20, 20, 1013.25, 35, 140)
    assert abs(result.magnitude) < 1e-12


def test_headspace_salinity_method_lowers_concentration():
    args = ("CH4", CH4_HS[0], CH4_AIR[0], 20, 20, 1013.25, VOL_HS, VOL_TOTAL)
    fresh = pyg.headspace_concentration(*args, method="Wiesenburg", salt=0)
    saline = pyg.headspace_concentration(*args, method="Wiesenburg", salt=35)
    # Lower solubility in salt water -> less gas left in the water after shaking
    assert saline.magnitude < fresh.magnitude


def test_headspace_invalid_volumes():
    with pytest.raises(ValueError):
        pyg.headspace_concentration("CH4", 5, 2, 20, 20, 1013.25, 140, 35)
