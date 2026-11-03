"""Tests for the water vapour correction (vapour=True) and vapour_pressure."""

import numpy as np
import pytest

import pyghgaq as pyg
from pyghgaq.functions.units import Q_


def yamamoto_ch4(temp_c, salt):
    """CH4 solubility (mol L-1 atm-1) from the Bunsen coefficient of Yamamoto et al. 1976."""
    t = (temp_c + 273.15) / 100
    bunsen = np.exp(-67.1962 + 99.1624 / t + 27.9015 * np.log(t)
                    + salt * (-0.072909 + 0.041674 * t - 0.0064603 * t**2))
    return bunsen / 22.414


def test_vapour_pressure_fresh_water():
    # Saturation vapour pressure of pure water at 25 degC: 31.69 hPa
    p = pyg.vapour_pressure(25)
    assert p.units == Q_(1, "hPa").units
    np.testing.assert_allclose(p.magnitude, 31.69, rtol=2e-3)


def test_vapour_pressure_lower_over_seawater():
    ratio = (pyg.vapour_pressure(25, 35) / pyg.vapour_pressure(25, 0)).magnitude
    # Raoult's law: about 1.9 % lower at 35 PSU
    assert 0.975 < ratio < 0.985


@pytest.mark.parametrize("temp", [0.0, 10.0, 20.0])
@pytest.mark.parametrize("salt", [0.0, 35.0])
def test_ch4_wiesenburg_vapour_free_matches_yamamoto(temp, salt):
    # W&G 1979 is the equilibrium with moist air at 1 atm; without the vapour
    # it is the Bunsen solubility of Yamamoto et al. 1976 (the two fits differ
    # by up to 0.3 %)
    from pyghgaq.main import _vapour_free_henry

    hcp = _vapour_free_henry("CH4", Q_(temp, "degC"), "Wiesenburg", {},
                             {"salt": salt, "catm": 1.9})
    np.testing.assert_allclose(hcp.to("mol L^-1 atm^-1").magnitude,
                               yamamoto_ch4(temp, salt), rtol=5e-3)


def test_vapour_false_is_unchanged():
    args = ("CO2", 600, 410, 10, 12, 1010, 35, 140)
    default = pyg.headspace_concentration(*args, method="Wiesenburg", salt=30)
    explicit = pyg.headspace_concentration(*args, method="Wiesenburg", salt=30,
                                           vapour=False)
    np.testing.assert_allclose(explicit.magnitude, default.magnitude, rtol=1e-12)


@pytest.mark.parametrize("gas, c_air", [("CH4", 1.9), ("CO2", 410.0)])
@pytest.mark.parametrize("rh_air", [0.0, 0.8, 1.0])
def test_vapour_mass_balance_zero_when_no_gas_added(gas, c_air, rh_air):
    # Gas-free seawater: the air put in is split between the headspace and the
    # water, so the water concentration before shaking must come out as zero.
    from pyghgaq.main import _vapour_free_henry

    temp_eq, temp_air, press, salt = 10.0, 5.0, 1005.0, 32.0
    v_hs, v_w = 35e-6, 105e-6
    r = 8.314462618
    hcp = _vapour_free_henry(gas, Q_(temp_eq, "degC"), "Wiesenburg", {},
                             {"salt": salt, "catm": c_air}).to("mol m^-3 Pa^-1").magnitude
    p_hs = (press - pyg.vapour_pressure(temp_eq, salt).magnitude) * 100
    p_air = (press - rh_air * pyg.vapour_pressure(temp_air).magnitude) * 100
    n_air = c_air * 1e-6 * p_air * v_hs / (r * (temp_air + 273.15))
    # n_air = x_hs p_hs (V_hs / (R T_eq) + H V_w)
    c_hs = n_air / (p_hs * (v_hs / (r * (temp_eq + 273.15)) + hcp * v_w)) * 1e6
    result = pyg.headspace_concentration(
        gas, c_hs, c_air, temp_eq, temp_air, press, 35, 140,
        method="Wiesenburg", salt=salt, vapour=True, rh_air=rh_air,
    )
    assert abs(result.magnitude) < 1e-9


def test_vapour_co2_matches_dry_pressure():
    # For CO2 (Weiss 1974, a coefficient per partial pressure) with the air
    # saturated at the shaking temperature, vapour=True is the default
    # calculation at the dry pressure P - pH2O
    temp, salt, press = 8.0, 33.0, 1012.0
    dry = press - pyg.vapour_pressure(temp, salt).magnitude
    with_vapour = pyg.headspace_concentration(
        "CO2", 350, 415, temp, temp, press, 35, 140,
        method="Wiesenburg", salt=salt, vapour=True, rh_air=1,
    )
    # The headspace after shaking is over seawater, the air before shaking over
    # fresh water (pH2O 1.8 % higher), so the two terms take different dry
    # pressures: headspace and water term (c_air = 0) + air term (c_hs = 0)
    dry_air = press - pyg.vapour_pressure(temp).magnitude
    hs_and_water = pyg.headspace_concentration(
        "CO2", 350, 0, temp, temp, dry, 35, 140, method="Wiesenburg", salt=salt,
    )
    air = pyg.headspace_concentration(
        "CO2", 0, 415, temp, temp, dry_air, 35, 140, method="Wiesenburg", salt=salt,
    )
    at_dry = hs_and_water + air
    np.testing.assert_allclose(with_vapour.magnitude, at_dry.magnitude, rtol=1e-10)


def test_csat_vapour_co2():
    temp, salt, press = 10.0, 35.0, 1013.25
    plain = pyg.csat("CO2", press, 410, temp, method="Wiesenburg", salt=salt)
    moist = pyg.csat("CO2", press, 410, temp, method="Wiesenburg", salt=salt, vapour=True)
    factor = 1 - pyg.vapour_pressure(temp, salt).magnitude / press
    np.testing.assert_allclose(moist.magnitude, plain.magnitude * factor, rtol=1e-10)


def test_csat_vapour_ch4_wiesenburg_at_1_atm_unchanged():
    # W&G already describe moist air at 1 atm: the vapour option must not
    # count the vapour twice
    plain = pyg.csat("CH4", 1013.25, 1.9, 10, method="Wiesenburg", salt=35)
    moist = pyg.csat("CH4", 1013.25, 1.9, 10, method="Wiesenburg", salt=35, vapour=True)
    np.testing.assert_allclose(moist.magnitude, plain.magnitude, rtol=1e-10)


@pytest.mark.parametrize("gas, c_hs, c_air, d_hs, d_air", [
    ("CH4", 6.85, 1.85, -66.7, -50.9),
    ("CO2", 1506.7, 403.4, -14.0, -8.65),
])
def test_headspace_d13c_vapour_equals_dry_pressure(gas, c_hs, c_air, d_hs, d_air):
    # With a coefficient per partial pressure (Sanders), fresh water, the air
    # saturated (rh_air=1) and temp_air = temp_eq, both partial pressures use
    # P - pH2O(temp): vapour=True is the default calculation at that pressure
    temp, press = 15.0, 1012.0
    dry = press - pyg.vapour_pressure(temp).magnitude
    args = (gas, c_hs, c_air, d_hs, d_air, temp, temp)
    with_vapour = pyg.headspace_d13c(*args, press, 35, 140, vapour=True, rh_air=1)
    at_dry = pyg.headspace_d13c(*args, dry, 35, 140)
    np.testing.assert_allclose(with_vapour.magnitude, at_dry.magnitude, rtol=1e-10)
    # Scaling all partial pressures by the same factor leaves d13C unchanged;
    # dry air before shaking (rh_air=0) does not, so d13C must change
    plain = pyg.headspace_d13c(*args, press, 35, 140)
    np.testing.assert_allclose(with_vapour.magnitude, plain.magnitude, rtol=1e-10)
    dry_air = pyg.headspace_d13c(*args, press, 35, 140, vapour=True, rh_air=0)
    assert abs(dry_air.magnitude - with_vapour.magnitude) > 1e-3


def test_vapour_salt_unit_override_reaches_vapour_pressure():
    # salt given as a mass fraction through units must give the same result as
    # the equivalent salinity in PSU, for the Henry's coefficient and for pH2O
    psu = {"method": "Wiesenburg", "salt": 35, "vapour": True}
    frac = {"method": "Wiesenburg", "salt": 0.035, "vapour": True,
            "units": {"salt": "dimensionless"}}
    for gas, catm in (("CO2", 410), ("CH4", 1.9)):
        np.testing.assert_allclose(
            pyg.csat(gas, 1005, catm, 10, **frac).magnitude,
            pyg.csat(gas, 1005, catm, 10, **psu).magnitude, rtol=1e-10,
        )
        args = (gas, catm * 3, catm, 10, 12, 1005, 35, 140)
        np.testing.assert_allclose(
            pyg.headspace_concentration(*args, **frac).magnitude,
            pyg.headspace_concentration(*args, **psu).magnitude, rtol=1e-10,
        )
