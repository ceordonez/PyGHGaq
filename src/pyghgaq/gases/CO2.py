import numpy as np
from pint import Quantity
from pint.facets.plain import PlainQuantity

from pyghgaq.functions.read import read_constant
from pyghgaq.functions.type_utils import Numeric
from pyghgaq.functions.units import Q_
from pyghgaq.registry.registry import enforce_units, register_hcp, register_schmidt


@register_hcp("CO2", "Weisenburg")
@enforce_units(temp="K", salt="PSU")
def hcp_sal_co2(
    varname: str, temp: Numeric, salt: Numeric, units: dict[str, str] = {}
) -> PlainQuantity:

    constant = read_constant()
    A = constant[varname]["A"]
    B = constant[varname]["B"]

    if isinstance(salt, Quantity):
        salt = salt.magnitude
    if isinstance(temp, Quantity):
        temp = temp.magnitude

    hcp: PlainQuantity= Q_(
        np.exp(
            +A[0]
            + A[1] * 100 / temp
            + A[2] * np.log(temp / 100)
            + salt * (B[0] + B[1] * temp / 100 + B[2] * (temp / 100) ** 2)
        ),
        "mol m^-3 atm^-1",
    )
    return hcp.to("mol m^-3 Pa^-1")


@register_hcp("CO2", "Sanders")
@enforce_units(temp="K")
def hcp_sanders(
    varname: str,
    temp: Quantity,
    units: dict[str, str] = {},
) -> Quantity:
    constant = read_constant()
    hcp25 = Q_(constant[varname]["H_T25"], "mol m^-3 Pa^-1")
    dlnHcpd1_T = Q_(constant[varname]["dlnHdT"], "K")
    hcp_t = hcp25 * np.exp(dlnHcpd1_T * (1 / temp - 1 / Q_(298.15, "K")))
    return hcp_t


@register_schmidt("CO2")
@enforce_units(temp="K")
def sch_number(temp: Numeric) -> float | np.ndarray:
    constant = read_constant()
    const = constant["CO2"]["SCH"][::-1]
    if isinstance(temp, Quantity):
        return np.polyval(const, temp.magnitude)
    else:
        return np.polyval(const, temp)
