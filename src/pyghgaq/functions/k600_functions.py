import numpy as np
from pyghgaq.registry.registry import register_k600


@register_k600("VP2013")
def k600_VP2013(u10_ms: np.ndarray | float, area_km2: float) -> np.ndarray | float:
    """Calculates gas transfer coefficient k600 from Vachon and Prairie 2013

    Parameters
    ----------
    u10_ms  : Wind velocity at 10m  (ms-1)
    area_km2 : Lake area (km2)

    Return
    ------
    k600 : velocity transfer coefficient (cmh-1)
    """
    return 2.51 + 1.48 * u10_ms + 0.39 * u10_ms * np.log10(area_km2)


@register_k600("MA2010-NB")
def k600_MA2010_NB(u10_ms: np.ndarray | float) -> np.ndarray | float:
    """Calculates gas transfer coefficient from McIntyre et al. 2010 negative bouyancy

    Parameters
    ----------
    u10_ms  : Wind velocity at 10m  (ms-1)

    Return
    ------
    k600 : velocity transfer coefficient (cmh-1)
    """
    return 2 + 2.04 * u10_ms


@register_k600("MA2010-PB")
def k600_MA2010_PB(u10_ms: np.ndarray | float) -> np.ndarray | float:
    """Calculates gas transfer coefficient from McIntyre et al. 2010 positive bouyancy

    Parameters
    ----------
    u10_ms  : Wind velocity at 10m  (ms-1)

    Return
    ------
    k600 : velocity transfer coefficient (cmh-1)
    """
    return 1.74 * u10_ms - 0.15


@register_k600("MA2010-MB")
def k600_MA2010_MB(u10_ms: np.ndarray | float) -> np.ndarray | float:
    """Calculates gas transfer coefficient from McIntyre et al. 2010 mixed model

    Parameters
    ----------
    u10_ms  : Wind velocity at 10m  (ms-1)

    Return
    ------
    k600 : velocity transfer coefficient (cmh-1)
    """
    return 2.25 * u10_ms + 0.16


@register_k600("CC1998")
def k600_CC1998(u10_ms: np.ndarray | float) -> np.ndarray | float:
    """Calculates gas transfer coefficient from Cole and Caraco 1998.

    Parameters
    ----------
    u10_ms  : Wind velocity at 10m  (ms-1)

    Return
    ------
    k600 : velocity transfer coefficient (cmh-1)
    """
    return 2.07 + 0.215 * u10_ms**1.7


def kx_k600(
    varname: str,
    temp_c: np.ndarray | float,
    k: np.ndarray | float,
    u10_ms: np.ndarray | float,
    a: int = 1,
) -> np.ndarray | float:
    """Calculates gas transfer coefficient kgas from k600 or vicecersa

    Parameters:
    ----------
    varname : gas name
    temp_c : water temperature (degC)
    k : gas transfer coefficient k600 or kgas
    u10_ms : wind speed at 10 m in (ms-1)
    a: if a ==  1 calculates kgas from k600
       if a == -1 calcualtes k600 from kgas

    Return
    ------
    kx : velocity transfer coefficient (md-1)
    """

    from pyghgaq.functions.functions import schmidt_number

    # Prairie and del Giorgo 2013
    if isinstance(u10_ms, np.ndarray):
        n = np.ones(len(u10_ms)) * 1 / 2
        n = np.where(u10_ms > 3.7, n, 2 / 3.0)
    else:
        n = 2 / 3.0
        if u10_ms < 3.7:
            n = 1 / 2.0

    sch = schmidt_number(varname, temp_c)
    return k * (600 / sch) ** (n * a)
