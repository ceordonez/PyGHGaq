import numpy as np
from functions.functions import read_constant
from registry.henry_registry import register_hcp

@register_hcp("CO2", "Weisenburg")
def hcp_sal_co2(
    varname: str,
    temp_c: np.ndarray | float,
    salt_psu: np.ndarray | float,
) -> np.ndarray | float:

    temp_k = temp_c + 273.15
    constant = read_constant()
    A = constant[varname]["A"]
    B = constant[varname]["B"]

    c_molm3atm = (
        np.exp(
            +A[0]
            + A[1] * 100 / temp_k
            + A[2] * np.log(temp_k / 100)
            + salt_psu * (B[0] + B[1] * temp_k / 100 + B[2] * (temp_k / 100) ** 2)
        )
        * 1000
    )
    hcpsalt_molm3Pa = c_molm3atm / 101325.0
    return hcpsalt_molm3Pa


@register_hcp("CO2", "Sanders")
def hcp_sanders(
    varname: str,
    temp_c: np.ndarray | float,
) -> np.ndarray | float:
    constant = read_constant()
    hcp25 = constant[varname]["H_T25"]
    dlnHcpd1_T = constant[varname]["dlnHdT"]
    hcp_t = hcp25 * np.exp(dlnHcpd1_T * (1 / (temp_c + 273.15) - 1 / 298.15))
    return hcp_t
