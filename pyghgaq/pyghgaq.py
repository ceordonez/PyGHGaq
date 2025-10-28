import numpy as np
import pandas as pd


def csat(
    atmpress_hpa: np.ndarray, catm_ppm: np.ndarray, hcp_molm3pa: np.ndarray
) -> np.ndarray:
    """
    Return
    -------
    Concentration of saturation in mmolm3
    """
    return atmpress_hpa * catm_ppm * hcp_molm3pa * 1e-2


def atm_diff_flux(
    csat_mmolm3: np.ndarray, cw_mmolm3: np.ndarray, k_md: np.ndarray
) -> np.ndarray:
    """
    Returns
    -------
    Diffusive flux to/from the atmosphere in mmolm2d
    """
    return k_md * (cw_mmolm3 - csat_mmolm3)


def hcp(
    varname: str,
    temp_c: pd.Series,
    sal_psu: pd.Series,
    catm_ppm: pd.Series,
    constant: dict,
    coeff: int = 1,
) -> pd.Series:
    """Calculate Henry's Coefficient

    Parameters
    ----------
    varname : string
        Could be ch4 or co2
    temp_c : pandas.Series
        Temperature of equilibrium.
    sal_psu : pandas.Series
        Water salinity in PSU
    constant : dict
        Constant data (defined in utils/constant.yml).
    coeff : float
        0 Use Sanders 2015 - Only consider correction by temperature (Use this for freshwaters).
        1 Use Wiss 1974 to estimate Henry's Coefficient for CO2 and Wiesenburg & Guinasso 1979 for CH4, considering Salinity and Temperature coerrections (default).

    Returns
    -------
    hcp_t : float
            Henry's coefficient for CH4 or CO2 in mmolm-3Pa-1
    """

    if coeff == 0:
        if varname == "ch4":
            hcp25 = constant["H_CH4_T25"]
            dlnHcpd1_T = constant["dlnH_CH4_d1T"]
        elif varname == "co2":
            hcp25 = constant["H_CO2_T25"]
            dlnHcpd1_T = constant["dlnH_CO2_d1T"]
        hcp_t = hcp25 * np.exp(dlnHcpd1_T * (1 / (temp_c + 273.15) - 1 / 298.15))
        return hcp_t

    elif coeff == 1:
        hcp_t = hcp_sal(varname, sal_psu, temp_c, catm_ppm)
        return hcp_t
    else:
        ## CHANGE TO RAISE AN ERROR
        return None


def hcp_sal(
    var: str, sal_psu: pd.Series, temp_c: pd.Series, catm_ppm: pd.Series
) -> pd.Series:
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
    temp_k = temp_c + 273.15
    fx = catm_ppm * 1e-6
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
        # A = [-417.5053, 599.8626, 380.3636, -62.0764]
        A = [-415.2807, 596.8104, 379.2599, -62.0757]
        B = [-0.05916, 0.032174, -0.0048198]

        c_molm3 = (
            np.exp(
                +np.log(fx)
                + A[0]
                + A[1] * 100 / temp_k
                + A[2] * np.log(temp_k / 100)
                + A[3] * (temp_k / 100)
                + sal_psu * (B[0] + B[1] * temp_k / 100 + B[2] * (temp_k / 100) ** 2)
            )
            * 1000
            / 1e9
        )
        hcpsalt_molm3Pa = c_molm3 / (101325 * fx)
    return hcpsalt_molm3Pa

