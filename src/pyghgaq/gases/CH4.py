import numpy as np
from pint import Quantity

from pyghgaq.functions.read import read_constant
from pyghgaq.functions.type_utils import Numeric
from pyghgaq.functions.units import Q_
from pyghgaq.registry.registry import enforce_units, register_hcp, register_schmidt


@register_hcp("CH4", "Weisenburg")
@enforce_units(temp="K", salt="PSU", catm="ppm")
def hcp_sal_ch4(
    varname: str,
    temp: Numeric,
    salt: Numeric,
    catm: Numeric,
    units: dict[str, str],
) -> Numeric:

    constant = read_constant()
    A = constant[varname]["A"]
    B = constant[varname]["B"]

    if isinstance(salt, Quantity):
        salt = salt.magnitude

    if isinstance(catm, Quantity):
        fx = catm.to_base_units().magnitude
    else:
        fx = catm * 1e-6
    if isinstance(temp, Quantity):
        temp = temp.magnitude

    c_molm3 = Q_(
        np.exp(
            +np.log(fx)
            + A[0]
            + A[1] * 100 / temp
            + A[2] * np.log(temp / 100)
            + A[3] * (temp / 100)
            + salt * (B[0] + B[1] * temp / 100 + B[2] * (temp / 100) ** 2)
        ),
        "nmol L^-1",
    )
    hcpsalt = c_molm3.to("mol m^-3") / (Q_(1, "atm").to("Pa") * fx)
    return hcpsalt


@register_hcp("CH4", "Sanders")
@enforce_units(temp="K")
def hcp_sanders(
    varname: str,
    temp: np.ndarray | Quantity | float,
    units: dict[str, str],
) -> np.ndarray | Quantity | float:

    constant = read_constant()
    hcp25 = Q_(constant[varname]["H_T25"], "mol m^-3 Pa^-1")
    dlnHcpd1_T = Q_(constant[varname]["dlnHdT"], "K")
    hcp_t = hcp25 * np.exp(dlnHcpd1_T * (1 / temp - 1 / Q_(298.15, "K")))
    return hcp_t


@register_schmidt("CH4")
@enforce_units(temp="K")
def sch_number(temp: Numeric) -> np.ndarray | float:
    constant = read_constant()
    const = constant["CH4"]["SCH"][::-1]
    if isinstance(temp, Quantity):
        return np.polyval(const, temp.magnitude)
    else:
        return np.polyval(const, temp)
