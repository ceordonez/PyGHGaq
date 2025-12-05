import numpy as np

import pyghgaq as pyg
from pyghgaq.functions.units import DEFAULT_INPUT_UNITS

print(DEFAULT_INPUT_UNITS)
DEFAULT_INPUT_UNITS["WindSpeed"] = "m/s"
DEFAULT_INPUT_UNITS["kgas"] = "m/d"

temp = np.arange(20, 25, 0.5)

temp = np.arange(1, 10, 0.5)
u = np.arange(1, 10, 0.5)
cw = np.arange(1, 10, 0.5)
cgas = np.arange(1, 10, 0.5)

value = pyg.k600_to_kgas("CH4", 1, 1, 1)
kgas = pyg.kgas(1, 1, 0.1)
atmf = pyg.atm_diff_flux(0, 1, 2, units={"kgas": "m/s"})
atmfa = pyg.atm_diff_flux(0, 1, 2)
k600 = pyg.kgas_to_k600("CH4", kgas, temp, u)
csat = pyg.csat("CH4", 1013.25, 2, 25)
u10 = pyg.uz_to_u10(u, 2, units={"windspeed": "km/h"})
print(value, kgas, atmf, atmfa, k600, csat, u10)
