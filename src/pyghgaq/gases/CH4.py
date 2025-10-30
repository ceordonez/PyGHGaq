import numpy as np
from pyghgaq.functions.read import read_constant
from pyghgaq.registry.registry import register_hcp, register_schmidt


@register_hcp("CH4", "Weisenburg")
def hcp_sal_ch4(
    varname: str,
    temp_c: np.ndarray | float,
    salt_psu: np.ndarray | float,
    catm_ppm: np.ndarray | float,
) -> np.ndarray | float:

    fx = catm_ppm * 1e-6
    temp_k = temp_c + 273.15
    constant = read_constant()
    A = constant[varname]["A"]
    B = constant[varname]["B"]

    c_molm3 = (
        np.exp(
            +np.log(fx)
            + A[0]
            + A[1] * 100 / temp_k
            + A[2] * np.log(temp_k / 100)
            + A[3] * (temp_k / 100)
            + salt_psu * (B[0] + B[1] * temp_k / 100 + B[2] * (temp_k / 100) ** 2)
        )
        * 1000
        / 1e9
    )
    hcpsalt_molm3Pa = c_molm3 / (101325.0 * fx)
    return hcpsalt_molm3Pa


@register_hcp("CH4", "Sanders")
def hcp_sanders(
    varname: str,
    temp_c: np.ndarray | float,
) -> np.ndarray | float:
    constant = read_constant()
    hcp25 = constant[varname]["H_T25"]
    dlnHcpd1_T = constant[varname]["dlnHdT"]
    hcp_t = hcp25 * np.exp(dlnHcpd1_T * (1 / (temp_c + 273.15) - 1 / 298.15))
    return hcp_t

@register_schmidt("CH4")
def schmit_number(temp_c):
    constant = read_constant()
    const = constant["CO2"]["SCH"][::-1]
    return np.polyval(const, temp_c)
