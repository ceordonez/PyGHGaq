import sys

sys.path.append("../../")

from pyghgaq.main import *

if __name__ == "__main__":

    temp = np.random.random(30)
    u = np.random.random(30)
    cw = np.random.random(30)
    cgas = np.random.random(30)
    # hcpch4_a = henry_coefficient("CH4", temp_c=temp)
    # hcpch4_b = henry_coefficient(
    #     "CH4", "Weisenburg", temp_c=temp, catm_ppm=2, salt_psu=0
    # )
    # k600_ms = k600("VP2013", u10_ms=1, area_km2=1)
    k600_ms = 10
    kgas_ms = k600_to_kgas("CO2", temp, k600_ms, u, ["degC", "m/h", "m/h"])
    print(kgas_ms)
    # co2flux = atm_diff_flux(cgas, cw, kgas_ms)
