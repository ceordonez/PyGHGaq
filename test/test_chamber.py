"""Tests for the chamber flux functions."""

import numpy as np
import pandas as pd
import pytest

import pyghgaq as pyg
from pyghgaq.functions.units import Q_


def test_chamber_flux_hand_calculation():
    # 1 ppm/s in 16.76 L over 0.126 m^2 at 30 degC and 1014 hPa
    slope, vol, area, temp, press = 1.0, 16.76, 0.126, 30.0, 1014.0
    n_air = press * 100 * vol * 1e-3 / (8.314462618 * (temp + 273.15))  # mol
    expected = slope * 1e-6 * n_air / area * 86400 * 1e3  # mmol m^-2 d^-1
    flux = pyg.chamber_flux(slope, vol, area, temp, press)
    assert flux.units == Q_(1, "mmol m^-2 d^-1").units
    assert flux.magnitude == pytest.approx(expected, rel=1e-12)


def test_chamber_flux_matches_lgr_surfflux_formula():
    # Formula used by PyChamberFlux (formerly LGR_SurfFlux) before it called pyghgaq (R = 8.31451)
    slope, vol, area, temp, press = 2.5e-4, 16.76, 0.126, 30.0, 1014.0
    molar_volume = 8.31451 * (temp + 273.15) / (press * 100) * 1000  # L/mol
    old = slope * vol / (molar_volume * area) * 86.4
    flux = pyg.chamber_flux(slope, vol, area, temp, press)
    assert flux.magnitude == pytest.approx(old, rel=1e-4)


def test_chamber_flux_units_and_arrays():
    slopes = np.array([1e-4, 2e-4, -1e-4])
    flux = pyg.chamber_flux(slopes, 16.76, 0.126, 30, 1014)
    flux_ml = pyg.chamber_flux(
        slopes, 16760, 1260, 30, 1014, units={"volume": "mL", "area": "cm^2"}
    )
    flux_q = pyg.chamber_flux(
        Q_(slopes * 60, "ppm/min"), 16.76, 0.126, Q_(303.15, "K"), Q_(101400, "Pa")
    )
    np.testing.assert_allclose(flux.magnitude, flux_ml.magnitude)
    np.testing.assert_allclose(flux.magnitude, flux_q.magnitude)
    assert flux.magnitude[2] < 0


def _synthetic(n=200, dt=0.995, slope=2.5e-4, intercept=1.99, noise=0.0):
    rng = np.random.default_rng(0)
    seconds = np.arange(n) * dt
    conc = intercept + slope * seconds + noise * rng.standard_normal(n)
    times = pd.Timestamp("2025-01-11 17:34:12.188") + pd.to_timedelta(seconds, "s")
    return times, seconds, conc


def test_chamber_slope_uses_real_time():
    times, seconds, conc = _synthetic()
    fit = pyg.chamber_slope(times, conc)
    assert fit.slope.to("ppm/s").magnitude == pytest.approx(2.5e-4, rel=1e-9)
    assert fit.intercept.to("ppm").magnitude == pytest.approx(1.99, rel=1e-9)
    assert fit.rvalue == pytest.approx(1.0)
    assert fit.n == 200
    # fitting against the sample number would be ~0.5 % low
    per_sample = np.polyfit(np.arange(200), conc, 1)[0]
    assert per_sample / fit.slope.magnitude == pytest.approx(0.995)


@pytest.mark.parametrize(
    "as_time",
    [
        lambda t, s: pd.Series(t),
        lambda t, s: t.to_numpy(),
        lambda t, s: list(t.to_pydatetime()),
        lambda t, s: s,
        lambda t, s: Q_(s / 60, "min"),
    ],
)
def test_chamber_slope_time_inputs(as_time):
    times, seconds, conc = _synthetic(noise=0.002)
    ref = pyg.chamber_slope(seconds, conc)
    fit = pyg.chamber_slope(as_time(times, seconds), pd.Series(conc))
    assert fit.slope.magnitude == pytest.approx(ref.slope.magnitude, rel=1e-9)


def test_chamber_slope_matches_numpy_regression():
    _, seconds, conc = _synthetic(noise=0.002)
    fit = pyg.chamber_slope(seconds, conc)
    (slope, intercept), cov = np.polyfit(seconds, conc, 1, cov="unscaled")
    resid = conc - (slope * seconds + intercept)
    stderr = np.sqrt(np.sum(resid**2) / (len(conc) - 2) * cov[0, 0])
    assert fit.slope.magnitude == pytest.approx(slope, rel=1e-9)
    assert fit.intercept.magnitude == pytest.approx(intercept, rel=1e-9)
    assert fit.stderr.magnitude == pytest.approx(stderr, rel=1e-9)
    assert fit.rvalue == pytest.approx(np.corrcoef(seconds, conc)[0, 1])


def test_chamber_slope_conc_units():
    _, seconds, conc = _synthetic()
    fit = pyg.chamber_slope(seconds, conc * 1000, units={"conc": "ppb"})
    assert fit.slope.to("ppm/s").magnitude == pytest.approx(2.5e-4, rel=1e-9)


def test_chamber_slope_errors():
    with pytest.raises(ValueError):
        pyg.chamber_slope([0, 1], [1.0, 2.0])
    with pytest.raises(ValueError):
        pyg.chamber_slope([0, 1, 2], [1.0, 2.0])
    with pytest.raises(ValueError):
        pyg.chamber_slope([5, 5, 5], [1.0, 2.0, 3.0])


def test_chamber_slope_into_chamber_flux():
    times, _, conc = _synthetic()
    fit = pyg.chamber_slope(times, conc)
    flux = pyg.chamber_flux(fit.slope, 16.76, 0.126, 30, 1014)
    assert flux.magnitude == pytest.approx(
        pyg.chamber_flux(2.5e-4, 16.76, 0.126, 30, 1014).magnitude
    )
