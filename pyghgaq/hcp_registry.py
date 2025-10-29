
import numpy as np
from functools import wraps
from typing import Any, Callable

type Data = [[np.ndarray|float, np.ndarray|float, np.ndarray|float]]
type HcpFn = Callable[[Data], np.ndarray|float]

exporters: dict[str, HcpFn] = {}

def register_hcp(varname: str):
    def decorator(func: HcpFn):
        @wraps(func)
        def wrapper(*args: Any, **kwards: Any) -> Any:
            return func(*args, **kwards)

        exporters[varname] = wrapper
        return wrapper

    return decorator

def hcpf(temp, fx, salt, varname):
    exporter = exporters.get(varname)
    return exporter(temp, fx, salt)

@register_hcp("CO2")
def hcp_sal_co2(
    temp_k: np.ndarray | float,
    fx: np.ndarray | float,
    salt_psu: np.ndarray | float,
) -> np.ndarray | float:

    print('CO2')
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
    print(hcpsalt_molm3Pa)
    return 1, hcpsalt_molm3Pa

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
if __name__ == "__main__":

    data : Data = (1, 1, 1)
    a = hcpf(1, 1, 1, 'CO2')
    print('value', a)
        
    # print(hcpf(data, 'CO2'))
