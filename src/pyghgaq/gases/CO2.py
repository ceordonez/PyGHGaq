import numpy as np
from pint import Quantity

from pyghgaq.functions.read import read_constant
from pyghgaq.functions.type_utils import Numeric
from pyghgaq.functions.units import Q_
from pyghgaq.registry.registry import enforce_units, register_hcp, register_schmidt


@enforce_units(temp="K", salt="psu")
@register_hcp("CO2", "Weisenburg")
def hcp_sal_co2(
    varname: str,
    temp: Numeric,
    salt: Numeric,
) -> Quantity:

    constant = read_constant()
    A = constant[varname]["A"]
    B = constant[varname]["B"]

    c_molm3atm = (
        np.exp(
            +A[0]
            + A[1] * 100 / temp
            + A[2] * np.log(temp / 100)
            + salt * (B[0] + B[1] * temp / 100 + B[2] * (temp / 100) ** 2)
        )
        * 1000
    )
    hcpsalt_molm3Pa = c_molm3atm / 101325.0
    return hcpsalt_molm3Pa


@register_hcp("CO2", "Sanders")
@enforce_units(temp="K")
def hcp_sanders(
    varname: str,
    temp: Quantity | float,
) -> Quantity:
    constant = read_constant()
    hcp25 = Q_(constant[varname]["H_T25"], "mol m^-3 Pa^-1")
    dlnHcpd1_T = Q_(constant[varname]["dlnHdT"], "K")
    hcp_t = hcp25 * np.exp(dlnHcpd1_T * (1 / temp - 1 / Q_(298.15, "K")))
    return hcp_t


@register_schmidt("CO2")
def schmidt_number(temp: Numeric):
    constant = read_constant()
    const = constant["CO2"]["SCH"][::-1]
    return np.polyval(const, temp.magnitude)
