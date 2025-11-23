import sys

sys.path.append("../../")

from pyghgaq.main import *

if __name__ == "__main__":

    from pyghgaq.functions.units import DEFAULT_INPUT_UNITS

    DEFAULT_INPUT_UNITS["WindSpeed"] = "m/s"
    DEFAULT_INPUT_UNITS["kgas"] = "m/d"

    print(DEFAULT_INPUT_UNITS)
    temp = np.arange(20, 25, 0.5)

    temp = np.arange(1, 10, 0.5)
    u = np.arange(1, 10, 0.5)
    cw = np.arange(1, 10, 0.5)
    cgas = np.arange(1, 10, 0.5)

    # hcpch4_a = henry_coefficient("CH4", temp)
    # hcpch4_b = henry_coefficient(
    #     "CH4", temp, method="Wiesenburg", salt=34, catm=1.41, units={"catm": "ppm"}
    # )
    # hcpco2_a = henry_coefficient("CO2", temp)
    # hcpco2_b = henry_coefficient("CO2", temp, "Wiesenburg", salt=0)
    #
    # sat = csat(800, 2, 3e-5)
    # fgas = kWindSpeedgas(1, 1, 0)
    kgas = kgas(1, 1, 0.1)
    atmf = atm_diff_flux(0, 1, 2, units={"kgas": "m/s"})
    atmfa = atm_diff_flux(0, 1, 2)
    total = atmf.to(atmfa.units) + atmfa
    k600 = kgas_to_k600('CH4', kgas, temp, u)
    csat = csat('CH4', 1013.25, 2, 25)
    u10 = uz_to_u10(u, 2)
    # fk600b = k600(u, "CC1998", units={"u10": "m/s"})
    # fk600 = k600(u, "CC1998")
    print(u10)
    print(k600)
    print(atmfa)
    print(csat)
    # print(hcpch4_a, hcpch4_b)
    # print(hcpco2_a, hcpco2_b)

    __import__("pdb").set_trace()
