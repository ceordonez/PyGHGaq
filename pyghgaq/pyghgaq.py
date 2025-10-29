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
    Concentration of saturation in mmolm3
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


def henry_coeff(
    varname: str,
    temp_c: np.ndarray | float,
    salt_inc: bool = False,
    sal_psu: np.ndarray | float = 0.0,
    catm_ppm: np.ndarray | float = 0.0,
) -> np.ndarray | float:
    """Calculate Henry's Coefficient

    Parameters
    ----------
    varname : Could be ch4 or co2
    temp_c : Temperature of equilibrium in °C.
    constant : Constant data (defined in utils/constant.yml).
    salt_inc : bool
    sal_psu : Water salinity in PSU

    Returns
    -------
    hcp_t : Henry's coefficient for CH4 or CO2 in mmolm-3Pa-1
    """

    constant = read_constant()
    if varname.upper() not in constant["GASES"]:
        raise ValueError(f"{varname} is not on the permited list {constant['GASES']}")

    if not salt_inc:
        # Use Sanders 2015 - Only consider correction by temperature (Use this for freshwaters).
        hcp_t = hcp_sanders(varname, temp_c)
    else:
        # Use Wiss 1974 to estimate Henry's Coefficient for CO2 and Wiesenburg & Guinasso 1979 for CH4, considering Salinity and Temperature coerrections.
        hcp_t = hcp_sal(varname, sal_psu, temp_c, catm_ppm)
    return hcp_t


def hcp_sanders(varname: str, temp_c: np.ndarray | float) -> np.ndarray | float:
    constant = read_constant()
    hcp25 = constant[varname]["H_T25"]
    dlnHcpd1_T = constant[varname]["dlnHdT"]
    hcp_t = hcp25 * np.exp(dlnHcpd1_T * (1 / (temp_c + 273.15) - 1 / 298.15))
    return hcp_t


def hcp_sal(
    var: str, sal_psu: np.ndarray | float, temp_c: np.ndarray | float, catm_ppm: np.ndarray | float
) -> np.ndarray | float:
    """Correct Henrys coefficient for CO2 (Weiss 1974) and CH4 (Wiesenburg and Guinasso 1979) considering salnity and temperature corrections.

    Parameters
    ----------
    var : string.
        Variable name (ch4, co2).
    sal_psu : pandas.Series or float.
        Salinity in PSU.
    temp_c : pandas.Series or float.
        Water temperature in deg C.

    Returns
    -------
    hcpsalt_mmolm3Pa : list
        Henry's coefficient estimated considering salinity and temperature corrections [molm-3Pa-1]
    """
    import henry_coeff

    temp_k = temp_c + 273.15
    fx = catm_ppm * 1e-6
    constant = read_constant()

    if var == "co2":  # folowing Weiss 1974
        A = [-58.0931, 90.5069, 22.2940]
        B = [0.027766, -0.025888, 0.0050578]
        hcpsalt_molm3Pa = (
            np.exp(
                +A[0]
                + A[1] * 100 / temp_k
                + A[2] * np.log(temp_k / 100)
                + sal_psu * (B[0] + B[1] * temp_k / 100 + B[2] * (temp_k / 100) ** 2)
            )
            * 1000
            / 101325
        )
    elif var == "ch4":  # folowing Wiesenburg and Guinasso 1979


if __name__ == "__main__":
    print(henry_coeff("CH4", np.array(25.0)))
