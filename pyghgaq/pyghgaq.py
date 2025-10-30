import importlib
import inspect
import pkgutil

import numpy as np


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
    return atmpress_hpa * catm_ppm * hcp_molm3pa * 1e-2


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


def henry_coefficient(
    varname: str, method: str = "Sanders", **kwards
) -> np.ndarray | float:

    import gases
    from registry.registry import exportershcp

    gasnames = []
    for _, modulename, _ in pkgutil.iter_modules(gases.__path__):
        gasnames.append(modulename)

    if varname not in gasnames:
        raise ValueError(
            f"{varname} is not included as gas to be analyze \b Valid gases are {gasnames}"
        )

    importlib.import_module(f"gases.{varname}")
    exporter = exportershcp.get(varname, {}).get(method)

    if exporter is None:
        raise ValueError(
            f"No method ='{method}' found for varname='{varname}' or '{varname}' is not included in functions \n Supported method for '{varname}' are {list(exportershcp.get(varname, {}).keys())}"
        )

    sig = inspect.signature(exporter)
    valid_kwards = {k: v for k, v in kwards.items() if k in sig.parameters}
    return exporter(varname, **valid_kwards)


def k600(k600_method: str = "MA2010-NB", **kwards):
    """Calculates gas transfer coefficient k600

    Return
    ------
    k600 : velocity transfer coefficient (md-1)
    """

    from registry.registry import exportersk600

    importlib.import_module("functions.k600_functions")

    exporter = exportersk600.get(k600_method)
    if exporter is None:
        raise ValueError(
            f"'{k600_method}' is not included as valid method in k600_models.py \n Supported methods are {list(exportersk600.keys())}"
        )

    sig = inspect.signature(exporter)
    valid_kwards = {k: v for k, v in kwards.items() if k in sig.parameters}
    return exporter(**valid_kwards) * 24 / 100


def k600_to_kx(
    varname: str,
    k600: np.ndarray | float,
    temp_c: np.ndarray | float,
    u10_ms: np.ndarray | float,
) -> np.ndarray | float:
    """Calculates gas transfer coefficient kgas from k600

    Parameters:
    ----------
    varname : gas name
    k600 : normalized gas transfer coefficient k600
    temp_c : water temperature (degC)
    u10_ms : wind speed at 10 m in (ms-1)

    Return
    ------
    kx : gas transfer coefficient
    """
    from functions.k600_functions import kx_k600

    return kx_k600(varname, temp_c, k600, u10_ms, 1)


def kx_to_k600(
    varname: str,
    kx: np.ndarray | float,
    temp_c: np.ndarray | float,
    u10_ms: np.ndarray | float,
) -> np.ndarray | float:
    """Calculates normalized gas transfer coefficient k600 from kgas

    Parameters:
    ----------
    varname : gas name
    kx : gas transfer coefficient for specific gas
    temp_c : water temperature (degC)
    u10_ms : wind speed at 10 m in (ms-1)

    Return
    ------
    k600 : normalized gas transfer coefficient
    """
    from functions.k600_functions import kx_k600

    return kx_k600(varname, temp_c, kx, u10_ms, -1)


def schmidt_number(varname: str, temp_c: np.ndarray | float):
    from functions.functions import schmidt_number

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
    hcpch4_a = henry_coefficient("CH4", temp_c=temp)
    hcpch4_b = henry_coefficient(
        "CH4", "Weisenburg", temp_c=temp, catm_ppm=2, salt_psu=0
    )
    k600_ms = k600("VP2013", u10_ms=1, area_km2=1)
    kgas_ms = k600_to_kx("CO2", temp, k600_ms, u)
    co2flux = atm_diff_flux(cgas, cw, kgas_ms)
    print(co2flux)
