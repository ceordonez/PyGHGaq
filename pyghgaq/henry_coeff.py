import numpy as np
from functions import *
from hcp_registry import *


@register_hcp("CH4")
def hcp_sal_ch4(
    temp_k: np.ndarray | float,
    fx: np.ndarray | float,
    salt_psu: np.ndarray | float,
) -> np.ndarray | float:

    # constant = read_constant()
    # # A = [-417.5053, 599.8626, 380.3636, -62.0764]
    # A = constant[varname]["A"]
    # B = constant[varname]["B"]
    #
    # c_molm3 = (
    #     np.exp(
    #         +np.log(fx)
    #         + A[0]
    #         + A[1] * 100 / temp_k
    #         + A[2] * np.log(temp_k / 100)
    #         + A[3] * (temp_k / 100)
    #         + salt_psu * (B[0] + B[1] * temp_k / 100 + B[2] * (temp_k / 100) ** 2)
    #     )
    #     * 1000
    #     / 1e9
    # )
    # hcpsalt_molm3Pa = c_molm3 / (101325.0 * fx)
    hcpsalt_molm3Pa = 2.0
    return hcpsalt_molm3Pa

@register_hcp("CO2")
def hcp_sal_co2(
    varname: str,
    temp_k: np.ndarray | float,
    fx: np.ndarray | float,
    salt_psu: np.ndarray | float,
) -> np.ndarray | float:

    # constant = read_constant()
    # # A = [-417.5053, 599.8626, 380.3636, -62.0764]
    # A = constant[varname]["A"]
    # B = constant[varname]["B"]
    #
    # c_molm3 = (
    #     np.exp(
    #         +np.log(fx)
    #         + A[0]
    #         + A[1] * 100 / temp_k
    #         + A[2] * np.log(temp_k / 100)
    #         + A[3] * (temp_k / 100)
    #         + salt_psu * (B[0] + B[1] * temp_k / 100 + B[2] * (temp_k / 100) ** 2)
    #     )
    #     * 1000
    #     / 1e9
    # )
    # hcpsalt_molm3Pa = c_molm3 / (101325.0 * fx)
    hcpsalt_molm3Pa = 1.0
    return hcpsalt_molm3Pa
