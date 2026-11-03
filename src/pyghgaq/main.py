import importlib
import inspect
import pkgutil
from typing import Any, NamedTuple

import numpy as np
from pint.facets.plain import PlainQuantity

from pyghgaq.functions.read import read_constant
from pyghgaq.functions.type_utils import Numeric
from pyghgaq.functions.units import Q_, get_default_units
from pyghgaq.registry.registry import _as_quantity, enforce_units


@enforce_units(
    unit_getter=get_default_units,
    unit_map={"windspeed": "WindSpeed", "height": "height"},
)
def uz_to_u10(
    windspeed: Numeric, height: float | int, units: dict[str, str] = {}
) -> PlainQuantity:
    constant = read_constant()
    cd = constant["Cd"]
    k = constant["k"]
    u10 = windspeed.to("m/s") * (1 + cd**0.5 / k * np.log(Q_(10, "m") / height))
    return u10


@enforce_units(
    unit_getter=get_default_units,
    unit_map={
        "atmpress": "Atm_Pressure",
        "catm": "Gas_Concentration",
        "temp": "Temperature",
    },
)
def csat(
    var: str,
    atmpress: Numeric,
    catm: Numeric,
    temp: Numeric,
    method: "str" = "Sanders",
    vapour: bool = False,
    units: dict[str, str] = {},
    **kwargs,
) -> PlainQuantity:
    """
    Concentration in equilibrium with air of mole fraction catm.

    With vapour=True, catm is a dry mole fraction and the air at the water
    surface is saturated with water vapour, so the partial pressure is
    catm (atmpress - pH2O(temp, salt)) (Weiss and Price 1980). The CH4
    "Wiesenburg" coefficient already gives the equilibrium with moist air at
    1 atm; without vapour it is used as is.

    Return
    -------
    Concentration of saturation in mmolm3
    """
    # catm is a named argument here, so forward it explicitly: methods such as
    # CH4 "Wiesenburg" require it (henry_coefficient drops it when unused).
    if vapour:
        hcp = _vapour_free_henry(var, temp, method, units, {"catm": catm, **kwargs})
        pressure = atmpress - _salt_vapour_pressure(temp, units, kwargs)
    else:
        hcp = henry_coefficient(var, temp, method, units, catm=catm, **kwargs)
        pressure = atmpress
    return pressure.to("Pa") * catm.to_base_units() * hcp.to("mmol m^-3 Pa^-1")


@enforce_units(
    unit_getter=get_default_units,
    unit_map={"temp": "Temperature", "salt": "Salinity"},
)
def vapour_pressure(
    temp: Numeric, salt: Numeric = 0, units: dict[str, str] = {}
) -> PlainQuantity:
    """Water vapour pressure at saturation over (sea)water, Weiss and Price 1980.

    Parameters
    ----------
    temp : water temperature (default degC)
    salt : salinity (default PSU), 0 for fresh water

    Return
    ------
    Vapour pressure (hPa)
    """
    c = read_constant()["WP80"]
    t = temp.to("K").magnitude / 100
    s = salt.to("PSU").magnitude
    return Q_(np.exp(c[0] + c[1] / t + c[2] * np.log(t) + c[3] * s), "atm").to("hPa")


def _salt_vapour_pressure(
    temp: PlainQuantity, units: dict[str, str], kwargs: dict
) -> PlainQuantity:
    """vapour_pressure for the salt in kwargs, keeping the caller's salt unit."""
    salt_units = {"salt": units["salt"]} if "salt" in units else {}
    return vapour_pressure(temp, kwargs.get("salt", 0), units=salt_units)


# Henry's coefficient methods that give the equilibrium with water-saturated
# air at 1 atm rather than with the dry gas (Wiesenburg and Guinasso 1979)
_MOIST_AIR_METHODS = {("CH4", "Wiesenburg")}


def _vapour_free_henry(
    varname: str, temp: PlainQuantity, method: str, units: dict[str, str], kwargs: dict
) -> PlainQuantity:
    """Henry's coefficient per unit partial pressure of the gas itself."""
    hcp = henry_coefficient(varname, temp, method, units, **kwargs)
    if (varname, method) in _MOIST_AIR_METHODS:
        ph2o = _salt_vapour_pressure(temp, units, kwargs)
        hcp = hcp / (1 - (ph2o / Q_(1, "atm")).to("dimensionless"))
    return hcp


@enforce_units(
    unit_getter=get_default_units,
    unit_map={
        "flux": "Diff_Flux",
        "cw": "Diss_Gas_Concentration",
        "csat": "Diss_Gas_Concentration",
    },
)
def kgas(
    flux: Numeric,
    cw: Numeric,
    csat: Numeric,
    units: dict[str, str] = {},
) -> PlainQuantity:
    """
    Return
    ------
    kgas from flux and concentrations measurements.
    """

    return flux.to("mmol/m^2/d") / (cw.to("mmol/m^3") - csat.to("mmol/m^3"))


@enforce_units(
    unit_getter=get_default_units,
    unit_map={"csat": "Diss_Gas_Concentration", "cw": "Diss_Gas_Concentration"},
)
def atm_diff_flux(
    csat: Numeric,
    cw: Numeric,
    kgas: Numeric,
    units: dict[str, str] = {},
) -> PlainQuantity:
    """
    Returns
    -------
    Diffusive flux to/from the atmosphere
    """
    return kgas.to("m/d") * (cw.to("mmol m^-3") - csat.to("mmol m^-3"))


@enforce_units(
    unit_getter=get_default_units,
    unit_map={"temp": "Temperature", "catm": "Gas_Concentration", "salt": "Salinity"},
)
def henry_coefficient(
    varname: str,
    temp: Numeric,
    method: str = "Sanders",
    units: dict[str, str] = {},
    **kwargs,
) -> PlainQuantity:

    from pyghgaq import gases
    from pyghgaq.registry.registry import exportershcp

    gasnames = []
    for _, modulename, _ in pkgutil.iter_modules(gases.__path__):
        gasnames.append(modulename)

    if varname not in gasnames:
        raise ValueError(
            f"{varname} is not included as gas to be analyze \b Valid gases are {gasnames}"
        )

    importlib.import_module(f"pyghgaq.gases.{varname}")
    exporter = exportershcp.get(varname, {}).get(method)

    if exporter is None:
        raise ValueError(
            f"No method ='{method}' found for varname='{varname}' or '{varname}' is not included in functions \n Supported method for '{varname}' are {list(exportershcp.get(varname, {}).keys())}"
        )
    sig = inspect.signature(exporter)
    valid_kwargs = {k: v for k, v in kwargs.items() if k in sig.parameters}
    valid_kwargs.update({"units": units})
    return exporter(varname, temp, **valid_kwargs)


@enforce_units(unit_getter=get_default_units, unit_map={"u10": "WindSpeed"})
def k600(
    u10: Numeric, k600_method: str = "MA2010-NB", units: dict[str, str] = {}, **kwargs
) -> PlainQuantity:
    """Calculates gas transfer coefficient k600

    Return
    ------
    k600 : velocity transfer coefficient (md-1)
    """

    from pyghgaq.registry.registry import exportersk600

    importlib.import_module("pyghgaq.functions.k600_functions")

    exporter = exportersk600.get(k600_method)
    if exporter is None:
        raise ValueError(
            f"'{k600_method}' is not included as valid method in k600_models.py \n Supported methods are {list(exportersk600.keys())}"
        )

    sig = inspect.signature(exporter)
    valid_kwargs = {k: v for k, v in kwargs.items() if k in sig.parameters}
    if "units" in sig.parameters:
        valid_kwargs["units"] = units
    return exporter(u10, **valid_kwargs).to("m/d")


@enforce_units(
    unit_getter=get_default_units, unit_map={"temp": "Temperature", "u10": "WindSpeed"}
)
def k600_to_kgas(
    varname: str,
    temp: Numeric,
    k: Numeric,
    u10: Numeric,
    units: dict[str, str] = {},
) -> PlainQuantity:
    """Calculates gas transfer coefficient kgas from k600

    Parameters:
    ----------
    varname : gas name
    k600 : normalized gas transfer coefficient k600
    temp : water temperature (degC)
    u10s : wind speed at 10 m in (ms-1)

    Return
    ------
    kgas : gas transfer coefficient
    """
    from pyghgaq.functions.k600_functions import kgas_k600

    return kgas_k600(varname, temp, k, u10, units, 1)


@enforce_units(
    unit_getter=get_default_units, unit_map={"u10": "WindSpeed", "temp": "Temperature"}
)
def kgas_to_k600(
    varname: str,
    kgas: Numeric,
    temp: Numeric,
    u10: Numeric,
    units: dict[str, str] = {},
) -> PlainQuantity:
    """Calculates normalized gas transfer coefficient k600 from kgas

    Parameters:
    ----------
    varname : gas name
    kgas : gas transfer coefficient for specific gas
    temp_c : water temperature (degC)
    u10_ms : wind speed at 10 m in (ms-1)

    Return
    ------
    k600 : normalized gas transfer coefficient
    """
    from pyghgaq.functions.k600_functions import kgas_k600

    return kgas_k600(varname, temp, kgas, u10, units, -1)


@enforce_units(unit_getter=get_default_units, unit_map={"temp": "Temperature"})
def schmidt_number(varname: str, temp: Numeric) -> float | np.ndarray:
    """Calculates Schmidt number for gases

    Parameters:
    ----------
    varname : gas name, e.g. "CO2", "CH4" (freshwater) or "CO2-salt",
              "CH4-salt" (seawater)
    temp : water temperature (default unit degC)

    Return
    ------
    Schmidt number
    """
    from pyghgaq.functions.functions import schmidt_number as _schmidt_number

    return _schmidt_number(varname, temp.to("degC"))


_HEADSPACE_UNIT_MAP = {
    "c_hs": "Gas_Concentration",
    "c_air": "Gas_Concentration",
    "temp_eq": "Temperature",
    "temp_air": "Temperature",
    "atmpress": "Atm_Pressure",
    "vol_hs": "Volume",
    "vol_total": "Volume",
    "rh_air": "Relative_Humidity",
}


def _headspace_moles(
    varname: str,
    c_hs: PlainQuantity,
    c_air: PlainQuantity,
    temp_eq: PlainQuantity,
    temp_air: PlainQuantity,
    atmpress: PlainQuantity,
    vol_hs: PlainQuantity,
    vol_total: PlainQuantity,
    method: str,
    vapour: bool,
    rh_air: PlainQuantity,
    units: dict[str, str],
    kwargs: dict,
) -> tuple[PlainQuantity, PlainQuantity, PlainQuantity, PlainQuantity]:
    """Mole balance of a headspace equilibration (inputs already Quantities).

    Returns
    -------
    n_hs : moles in the headspace after shaking
    n_air : moles in the headspace before shaking
    n_aq : moles in the water after shaking
    vol_water : water volume
    """
    vol_water = vol_total.to("m^3") - vol_hs.to("m^3")
    if np.any(vol_water.magnitude <= 0):
        raise ValueError("vol_total must be larger than vol_hs")

    r_gas = Q_(read_constant()["R"], "J mol^-1 K^-1")
    # catm only matters for methods that need it (CH4 "Wiesenburg"); the
    # resulting Henry's coefficient does not depend on its value.
    hcp_kwargs = {"catm": c_hs, **kwargs}
    if vapour:
        # Dry mole fractions: the headspace is saturated with water vapour at
        # temp_eq after shaking; the air before shaking has rh_air at temp_air
        p_dry_hs = atmpress - _salt_vapour_pressure(temp_eq, units, kwargs)
        p_dry_air = atmpress - rh_air.to("dimensionless") * vapour_pressure(temp_air)
        hcp = _vapour_free_henry(varname, temp_eq, method, units, hcp_kwargs)
    else:
        p_dry_hs = p_dry_air = atmpress
        hcp = henry_coefficient(varname, temp_eq, method, units, **hcp_kwargs)
    p_hs = c_hs.to("dimensionless") * p_dry_hs.to("Pa")
    p_air = c_air.to("dimensionless") * p_dry_air.to("Pa")
    n_hs = (p_hs * vol_hs / (r_gas * temp_eq.to("K"))).to("mol")
    n_air = (p_air * vol_hs / (r_gas * temp_air.to("K"))).to("mol")
    n_aq = (hcp * p_hs * vol_water).to("mol")
    return n_hs, n_air, n_aq, vol_water


@enforce_units(unit_getter=get_default_units, unit_map=_HEADSPACE_UNIT_MAP)
def headspace_concentration(
    varname: str,
    c_hs: Numeric,
    c_air: Numeric,
    temp_eq: Numeric,
    temp_air: Numeric,
    atmpress: Numeric,
    vol_hs: Numeric,
    vol_total: Numeric,
    method: str = "Sanders",
    vapour: bool = False,
    rh_air: Numeric = 1.0,
    units: dict[str, str] = {},
    **kwargs,
) -> PlainQuantity:
    """Dissolved gas concentration from a headspace equilibration.

    The water sample (vol_total - vol_hs) is shaken with a headspace of
    volume vol_hs filled with air (or another gas) of known concentration.
    The concentration in the water before shaking is obtained from the mole
    balance:

        c_w = (n_hs + n_aq - n_air) / V_water

    where n_hs = p_hs V_hs / (R T_eq) is the gas in the headspace after
    shaking, n_aq = H_cp p_hs V_water the gas left in the water after
    shaking, and n_air = p_air V_hs / (R T_air) the gas in the headspace
    before shaking. Partial pressures are mole fraction x atmpress, or, with
    vapour=True, dry mole fraction x (atmpress - pH2O): pH2O at saturation
    at temp_eq (and salt) after shaking, and rh_air x pH2O at temp_air (fresh
    water) before shaking (Weiss and Price 1980). With vapour=True the CH4
    "Wiesenburg" coefficient, which includes the water vapour, is converted to
    the coefficient per partial pressure of CH4.

    Parameters
    ----------
    varname : gas name ("CH4" or "CO2")
    c_hs : headspace concentration after shaking (default ppm)
    c_air : headspace concentration before shaking (default ppm)
    temp_eq : equilibration (water) temperature (default degC)
    temp_air : air temperature when the headspace was filled (default degC)
    atmpress : atmospheric pressure (default hPa)
    vol_hs : headspace volume (default mL)
    vol_total : total volume of the bottle/syringe (default mL)
    method : Henry's coefficient method, see henry_coefficient
    vapour : correct the partial pressures for water vapour (c_hs and c_air
             are then dry mole fractions, as from a Picarro)
    rh_air : relative humidity (0-1) of the air before shaking, with vapour
    **kwargs : extra arguments for the Henry's coefficient method,
               e.g. salt for "Wiesenburg" (also used for pH2O)

    Return
    ------
    Dissolved gas concentration before shaking (mmol m^-3)
    """
    n_hs, n_air, n_aq, vol_water = _headspace_moles(
        varname, c_hs, c_air, temp_eq, temp_air, atmpress, vol_hs, vol_total,
        method, vapour, rh_air, units, kwargs,
    )
    return ((n_hs + n_aq - n_air) / vol_water).to("mmol m^-3")


@enforce_units(
    unit_getter=get_default_units,
    unit_map={
        **_HEADSPACE_UNIT_MAP,
        "d13c_hs": "Isotope_Signature",
        "d13c_air": "Isotope_Signature",
    },
)
def headspace_d13c(
    varname: str,
    c_hs: Numeric,
    c_air: Numeric,
    d13c_hs: Numeric,
    d13c_air: Numeric,
    temp_eq: Numeric,
    temp_air: Numeric,
    atmpress: Numeric,
    vol_hs: Numeric,
    vol_total: Numeric,
    method: str = "Sanders",
    vapour: bool = False,
    rh_air: Numeric = 1.0,
    units: dict[str, str] = {},
    **kwargs,
) -> PlainQuantity:
    """Isotopic signature (d13C) of the dissolved gas from a headspace equilibration.

    13C and 12C mass balance equivalent to headspace_concentration. Assumes
    no isotopic fractionation between the water and the headspace after
    shaking, i.e. the gas left in the water has the headspace signature
    d13c_hs. Signatures are relative to VPDB (constant RVDPD).

    Parameters
    ----------
    varname : gas name ("CH4" or "CO2")
    c_hs, c_air : headspace concentration after / before shaking (default ppm)
    d13c_hs, d13c_air : d13C of the headspace after / before shaking
                        (default permille)
    temp_eq, temp_air, atmpress, vol_hs, vol_total, method, vapour, rh_air,
    **kwargs : see headspace_concentration

    Return
    ------
    d13C of the dissolved gas before shaking (permille)
    """
    n_hs, n_air, n_aq, _ = _headspace_moles(
        varname, c_hs, c_air, temp_eq, temp_air, atmpress, vol_hs, vol_total,
        method, vapour, rh_air, units, kwargs,
    )
    r_std = read_constant()["RVDPD"]

    def split(n, d13c):
        ratio = r_std * (d13c.to("permille").magnitude / 1000 + 1)
        n12 = n / (1 + ratio)
        return n - n12, n12

    hs13, hs12 = split(n_hs, d13c_hs)
    air13, air12 = split(n_air, d13c_air)
    aq13, aq12 = split(n_aq, d13c_hs)
    w13 = hs13 + aq13 - air13
    w12 = hs12 + aq12 - air12
    ratio_w = (w13 / w12).to("dimensionless").magnitude
    return Q_((ratio_w / r_std - 1) * 1000, "permille")


class SlopeFit(NamedTuple):
    """Least-squares line fitted by chamber_slope."""

    slope: PlainQuantity
    intercept: PlainQuantity
    rvalue: float
    stderr: PlainQuantity
    n: int


def _elapsed_seconds(time: Any) -> np.ndarray:
    """Seconds since the first sample, for datetimes or numbers (seconds)."""
    if isinstance(time, PlainQuantity):
        return time.to("s").magnitude.astype(float)
    t = np.asarray(time)
    if t.dtype == object:
        # e.g. lists of datetime.datetime or pandas Timestamps
        t = np.asarray(t, dtype="datetime64[ns]")
    if np.issubdtype(t.dtype, np.datetime64):
        return (t - t[0]) / np.timedelta64(1, "s")
    return t.astype(float)


def chamber_slope(
    time: Any, conc: Numeric, units: dict[str, str] = {}
) -> SlopeFit:
    """Linear rate of change of the gas concentration inside a chamber.

    The line is fitted against the real elapsed time, so samples that are
    not exactly 1 s apart (e.g. an LGR/Picarro logging every ~0.995 s, fast
    mode, or gaps in the record) give the correct slope per second.

    Parameters
    ----------
    time : sample times, as datetimes (datetime64, pandas Series or
           DatetimeIndex, datetime objects) or as numbers in seconds
    conc : gas concentration (default ppm, key "conc" in units)

    Return
    ------
    SlopeFit(slope, intercept, rvalue, stderr, n): slope and its standard
    error in ppm s^-1, intercept (concentration at the first sample) in ppm,
    Pearson correlation coefficient and number of samples.
    """
    x = _elapsed_seconds(time)
    unit = units.get("conc", get_default_units()["Gas_Concentration"])
    y = _as_quantity(conc, unit).to("ppm").magnitude
    if x.shape != y.shape:
        raise ValueError("time and conc must have the same length")
    n = x.size
    if n < 3:
        raise ValueError("chamber_slope needs at least 3 samples")

    dx = x - x.mean()
    dy = y - y.mean()
    sxx = np.sum(dx**2)
    syy = np.sum(dy**2)
    if sxx == 0:
        raise ValueError("all time values are equal")
    slope = np.sum(dx * dy) / sxx
    intercept = y.mean() - slope * x.mean()
    rvalue = float(np.sum(dx * dy) / np.sqrt(sxx * syy)) if syy > 0 else 0.0
    resid = dy - slope * dx
    stderr = np.sqrt(np.sum(resid**2) / (n - 2) / sxx)
    return SlopeFit(
        Q_(slope, "ppm/s"), Q_(intercept, "ppm"), rvalue, Q_(stderr, "ppm/s"), n
    )


@enforce_units(
    unit_getter=get_default_units,
    unit_map={
        "slope": "Conc_Slope",
        "volume": "Chamber_Volume",
        "area": "Chamber_Area",
        "temp": "Temperature",
        "atmpress": "Atm_Pressure",
    },
)
def chamber_flux(
    slope: Numeric,
    volume: Numeric,
    area: Numeric,
    temp: Numeric,
    atmpress: Numeric,
    units: dict[str, str] = {},
) -> PlainQuantity:
    """Flux from the rate of change of concentration in a (floating) chamber.

        F = (dC/dt) P V / (R T A)

    where P V / (R T) are the moles of air in the chamber (ideal gas).

    Parameters
    ----------
    slope : rate of change of the mole fraction in the chamber
            (default ppm s^-1), e.g. chamber_slope(...).slope
    volume : chamber volume, including tubing and analyser (default L)
    area : water surface covered by the chamber (default m^2)
    temp : air temperature in the chamber (default degC)
    atmpress : pressure in the chamber (default hPa)

    Return
    ------
    Flux (mmol m^-2 d^-1), positive = emission to the atmosphere
    """
    r_gas = Q_(read_constant()["R"], "J mol^-1 K^-1")
    n_air = atmpress.to("Pa") * volume.to("m^3") / (r_gas * temp.to("K"))
    return (slope.to("1/s") * n_air / area.to("m^2")).to("mmol m^-2 d^-1")


if __name__ == "__main__":

    temp = np.random.random(30)
    u = np.random.random(30)
    cw = np.random.random(30)
    cgas = np.random.random(30)
    # hcpch4_a = henry_coefficient("CH4", temp_c=temp)
    # hcpch4_b = henry_coefficient(
    #     "CH4", "Weisenburg", temp_c=temp, catm_ppm=2, salt_psu=0
    # )
    # co2flux = atm_diff_flux(cgas, cw, kgas_ms)
