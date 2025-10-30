import importlib
import inspect

import numpy as np
from functions import *


def csat(
    atmpress_hpa: np.ndarray | float,
    catm_ppm: np.ndarray | float,
    hcp_molm3pa: np.ndarray | float,
) -> np.ndarray | float:
    """
    Return
    -------
    Concentration of saturation
    """
    return atmpress_hpa * catm_ppm * hcp_molm3pa * 1e-2


def atm_diff_flux(
    csat_mmolm3: np.ndarray | float,
    cw_mmolm3: np.ndarray | float,
    k_md: np.ndarray | float,
) -> np.ndarray | float:
    """
    Returns
    -------
    Diffusive flux to/from the atmosphere in mmolm2d
    """
    return k_md * (cw_mmolm3 - csat_mmolm3)


def henry_coefficient(
    varname: str, method: str = "Sanders", **kwards
) -> np.ndarray | float:
    """
    Returns
    -------
    Henry's coefficient in molm3Pa
    """

    from registry.henry_registry import exporters

    importlib.import_module(f"gases.{varname}")
    exporter = exporters.get(varname, {}).get(method)
    if exporter is None:
        raise ValueError(
            f"No method ='{method}' found for varname='{varname}' or '{varname}' is not included in functions"
        )

    sig = inspect.signature(exporter)
    valid_kwards = {k: v for k, v in kwards.items() if k in sig.parameters}
    return exporter(varname, **valid_kwards)


if __name__ == "__main__":

    temp = np.arange(1, 30, 1)
    a = henry_coefficient("CH4", temp_c=temp)
    b = henry_coefficient("CH4", "Weisenburg", temp_c=temp, catm_ppm=2, salt_psu=0)
    print(a, b)
