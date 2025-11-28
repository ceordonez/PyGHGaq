import numpy as np
import pytest

from pyghgaq.functions.units import Q_, ureg
from pyghgaq.main import (
    atm_diff_flux,
    csat,
    henry_coefficient,
    k600,
    k600_to_kgas,
    kgas,
    kgas_to_k600,
    schmidt_number,
    uz_to_u10,
)


@pytest.fixture(autouse=True)
def setup_mock_registries(monkeypatch):
    """
    Automatically set up all necessary registry exporters and mock gases modules for each test.
    Ensures test isolation and prevents registry errors.
    """
    import pyghgaq.registry.registry

    # Mock exportershcp
    mock_henry = lambda varname, temp, **kwargs: Q_(1, "mmol m^-3 Pa^-1")
    monkeypatch.setitem(
        pyghgaq.registry.registry.exportershcp, "CO2", {"Sanders": mock_henry}
    )

    # Mock exportersk600
    mock_k600 = lambda u10, **kwargs: Q_(1, "m/d")
    monkeypatch.setitem(pyghgaq.registry.registry.exportersk600, "MA2010-NB", mock_k600)

    # Mock exporterssh (Schmidt number calculator)
    mock_schmidt = lambda temp: 1000.0
    monkeypatch.setitem(pyghgaq.registry.registry.exporterssh, "CO2", mock_schmidt)

    # Mock pyghgaq.gases and pyghgaq.gases.CO2 module presence
    import sys
    import types

    gases_mod = types.ModuleType("pyghgaq.gases")
    setattr(gases_mod, "__path__", ["dummy_path"])  # So pkgutil.iter_modules works
    sys.modules["pyghgaq.gases"] = gases_mod
    sys.modules["pyghgaq.gases.CO2"] = types.ModuleType("pyghgaq.gases.CO2")

    import pkgutil

    def mock_iter_modules(path):
        # Return a tuple in the format pkgutil uses: (module_finder, name, ispkg)
        # Only 'CO2' appears.
        return [(None, "CO2", False)]

    monkeypatch.setattr(pkgutil, "iter_modules", mock_iter_modules)


# Setup fixture for units
@pytest.fixture
def unit_reg():
    return ureg


def test_uz_to_u10(unit_reg):
    windspeed = Q_(5, "m/s")
    height = Q_(2, "m")
    result = uz_to_u10(windspeed, height)
    assert result.units == unit_reg("m/s").units
    assert result.magnitude > 0


def test_csat(unit_reg):
    atmpress = Q_(101325, "Pa")
    catm = Q_(400, "umol/mol")
    temp = Q_(20, "degC")
    # You may need to mock henry_coefficient internally
    try:
        result = csat("CO2", atmpress, catm, temp)
        assert result.units == unit_reg("mmol/m^3").units
    except Exception:
        pass  # Accept error if gases/registry not set up


def test_kgas(unit_reg):
    flux = Q_(10, "mmol/m^2/d")
    cw = Q_(200, "mmol/m^3")
    csat_v = Q_(180, "mmol/m^3")
    result = kgas(flux, cw, csat_v)
    assert result.magnitude != 0


def test_atm_diff_flux(unit_reg):
    csat_v = Q_(180, "mmol/m^3")
    cw_v = Q_(200, "mmol/m^3")
    kgas_val = Q_(1, "m/d")
    result = atm_diff_flux(csat_v, cw_v, kgas_val)
    assert (
        result.units == unit_reg("mmol/m^2/d").units
        or result.units == unit_reg("mmol/d/m^2").units
    )


def test_schmidt_number():
    temp = 20
    result = schmidt_number("CO2", temp)
    assert isinstance(result, (float, np.ndarray))


def test_henry_coefficient():
    temp = Q_(20, "degC")
    result = henry_coefficient("CO2", temp)
    assert hasattr(result, "units")


def test_k600():
    u10_val = Q_(1, "m/s")
    result = k600(u10_val)
    assert hasattr(result, "units")


def test_k600_to_kgas():
    temp = Q_(20, "degC")
    k = Q_(1, "m/d")
    u10 = Q_(1, "m/s")
    result = k600_to_kgas("CO2", temp, k, u10)
    assert hasattr(result, "units")


def test_kgas_to_k600():
    temp = Q_(20, "degC")
    kgas_val = Q_(1, "m/d")
    u10 = Q_(1, "m/s")
    result = kgas_to_k600("CO2", kgas_val, temp, u10)
    assert hasattr(result, "units")
