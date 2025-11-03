import importlib
import inspect
import pkgutil

import numpy as np

from pyghgaq.functions.type_utils import Numeric
from pyghgaq.registry.registry import warn_default_units


def csat(
    atmpress_hpa: np.ndarray | float,
    catm_ppm: np.ndarray | float,
    hcp_molm3pa: np.ndarray | float,
) -> np.ndarray | float:
    """
    Return
    -------
    Concentration of saturation in mmolm3
    """
    return atmpress_hpa * catm_ppm * hcp_molm3pa * 1e-1


def kgas(
    flux: np.ndarray | float,
    cw: np.ndarray | float,
    csat: np.ndarray | float,
) -> np.ndarray | float:
    """
    Return
    ------
    kgas from flux and concentrations measurements.
    """

    return flux / (cw - csat)


def atm_diff_flux(
    csat: np.ndarray | float,
    cw: np.ndarray | float,
    kgas: np.ndarray | float,
) -> np.ndarray | float:
    """
    Returns
    -------
    Diffusive flux to/from the atmosphere
    """
    return kgas * (cw - csat)


@warn_default_units
def henry_coefficient(
    varname: str,
    temp: np.ndarray | float,
    method: str = "Sanders",
    units: dict[str, str] = {"temp": "degC", "salt": "PSU", "catm": "ppm"},
    **kwargs,
) -> Numeric:

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


def k600(k600_method: str = "MA2010-NB", **kwargs):
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
    return exporter(**valid_kwargs) * 24 / 100


@warn_default_units
def k600_to_kgas(
    varname: str,
    temp: np.ndarray | float,
    k: np.ndarray | float,
    u10: np.ndarray | float,
    units: dict[str, str] = {"temp": "degC", "k600": "m/s", "u10": "m/s"},
) -> Numeric:
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

    units["k"] = units.pop("k600")
    return kgas_k600(varname, temp, k, u10, units, 1)


@warn_default_units
def kgas_to_k600(
    varname: str,
    kgas: np.ndarray | float,
    temp_c: np.ndarray | float,
    u10_ms: np.ndarray | float,
) -> np.ndarray | float:
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

    return kgas_k600(varname, temp_c, kgas, u10_ms, -1)


def schmidt_number(varname: str, temp_c: np.ndarray | float):
    from pyghgaq.functions.functions import schmidt_number

    """Calculates Schmidt number for gases

    Parameters:
    ----------
    varname : gas name
    temp_c : water temperature (degC)

    Return
    ------
    Schmidt number
    """

    return schmidt_number(varname, temp_c)


if __name__ == "__main__":

    temp = np.random.random(30)
    u = np.random.random(30)
    cw = np.random.random(30)
    cgas = np.random.random(30)
    # hcpch4_a = henry_coefficient("CH4", temp_c=temp)
    # hcpch4_b = henry_coefficient(
    #     "CH4", "Weisenburg", temp_c=temp, catm_ppm=2, salt_psu=0
    # )
    k600_ms = k600("VP2013", u10_ms=1, area_km2=1)
    kgas_ms = k600_to_kgas("CO2", temp, k600_ms, u, ["m/s", "degC", "m/s"])
    # co2flux = atm_diff_flux(cgas, cw, kgas_ms)
