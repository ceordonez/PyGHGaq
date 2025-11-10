import sys

sys.path.append("../../")

from pyghgaq.main import *

if __name__ == "__main__":

    temp = np.arange(20, 25, 0.5)

    temp = np.arange(1, 10, 0.5)
    u = np.arange(1, 10, 0.5)
    cw = np.arange(1, 10, 0.5)
    cgas = np.arange(1, 10, 0.5)

    hcpch4_a = henry_coefficient("CH4", temp)
    hcpch4_b = henry_coefficient("CH4", temp, method='Weisenburg', salt=34, catm=1.41, units={"catm": "ppm"})
    sat = csat(800, 2, 3E-5)
    fgas = kgas(1, 1, 0)
    atmf = atm_diff_flux(0, 1, 2)
    fk600 = k600(u, 'CC1998')


    print(hcpch4_a, hcpch4_b)
    hcpco2_a = henry_coefficient("CO2", temp)
    hcpco2_b = henry_coefficient("CO2", temp, 'Weisenburg', salt=1, catm=1)
    print(hcpco2_a, hcpco2_b)
    __import__('pdb').set_trace()
    # hcp_sal_ch4("CH4", 25, 0, 2, catm_units='ppb')
    # __import__("pdb").set_trace()
    # k600_ms = k600("VP2013", u10_ms=1, area_km2=1)
    # k600_ms = 10
    #
    # kgas_ms_a = k600_to_kgas(
    #     "CO2", temp, k600_ms, u
    # )
    # kgas_ms_b = k600_to_kgas(
    #     "CO2", temp, k600_ms, u, units={"temp": "degC", "k600": "m/h", "u10": "m/h"}
    # )
    # print(kgas_ms_a, kgas_ms_b)
    # co2flux = atm_diff_flux(cgas, cw, kgas_ms)
