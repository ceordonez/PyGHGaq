# PyGHGaq

This package contains methods to estimate GHG (CH<sub>4</sub> and CO<sub>2</sub>) emissions and concentrations in aquatic ecosystems.

Implemented methods:
* Diffusive flux through the water-air interface (`atm_diff_flux`) and gas transfer coefficient from measured fluxes (`kgas`).
* Gas transfer coefficient k<sub>600</sub> (`k600`) using different models:
    * MacIntyre et al. 2010: negative buoyancy (`"MA2010-NB"`, default), positive buoyancy (`"MA2010-PB"`) and all buoyancy conditions (`"MA2010-MB"`).
    * Cole and Caraco 1998 (`"CC1998"`).
    * Vachon and Prairie 2013 (`"VP2013"`, requires the lake `area`).
* Conversion from k<sub>600</sub> to k<sub>gas</sub> and vice versa (`k600_to_kgas`, `kgas_to_k600`).
* Schmidt number for freshwater (`"CO2"`, `"CH4"`) and seawater (`"CO2-salt"`, `"CH4-salt"`), Wanninkhof 1992 (`schmidt_number`).
* Henry's coefficient corrected by temperature (`henry_coefficient`):
    * `"Sanders"` (default): temperature correction only, Sander 2023.
    * `"Wiesenburg"`: temperature and salinity correction. Uses Wiesenburg and Guinasso 1979 for CH<sub>4</sub> (requires `salt` and `catm`) and Weiss 1974 for CO<sub>2</sub> (requires `salt`).
* Saturation concentration in water (`csat`).
* Dissolved gas concentration (`headspace_concentration`) and its δ<sup>13</sup>C (`headspace_d13c`) from headspace equilibrations.
* Wind speed at height z to wind speed at 10 m (`uz_to_u10`).

## Installation
```sh
pip install pyghgaq
```

## Usage

```python
import numpy as np
import pyghgaq as pyg

temp = np.array([10.0, 15.0, 20.0])  # water temperature (degC)
u10 = np.array([1.0, 3.0, 5.0])      # wind speed at 10 m (m/s)
cw = np.array([0.5, 0.8, 1.2])       # dissolved CH4 in water (mmol m^-3)

# Wind speed measured at 2 m -> 10 m
u10 = pyg.uz_to_u10(u10, 2)

# Saturation concentration for 1.9 ppm atmospheric CH4 at 1013.25 hPa
cs = pyg.csat("CH4", 1013.25, 1.9, temp)
cs_salt = pyg.csat("CH4", 1013.25, 1.9, temp, method="Wiesenburg", salt=0)

# Henry's coefficient
hcp = pyg.henry_coefficient("CH4", temp)

# Gas transfer coefficient
k600 = pyg.k600(u10)                            # MacIntyre et al. 2010 (negative buoyancy)
k600_vp = pyg.k600(u10, "VP2013", area=1.5)     # lake area in km^2
kch4 = pyg.k600_to_kgas("CH4", temp, k600, u10)

# Diffusive flux (positive = emission to the atmosphere)
flux = pyg.atm_diff_flux(cs, cw, kch4)
```

### Headspace method

A water sample is shaken with a headspace of known gas composition, and the headspace is measured after equilibration (for example with a Picarro). The concentration in the water before shaking follows from the mole balance

$$c_w = \frac{n_{HS} + n_{aq} - n_{air}}{V_w}, \qquad n_{HS} = \frac{p_{HS} V_{HS}}{R T_{eq}}, \qquad n_{aq} = H^{cp} p_{HS} V_w, \qquad n_{air} = \frac{p_{air} V_{HS}}{R T_{air}}$$

```python
import pyghgaq as pyg

# headspace after (6.85 ppm) and before (1.85 ppm) shaking, equilibration
# temperature 26.1 degC, air temperature 26.2 degC, 1014.7 hPa,
# 35 mL headspace in a 140 mL bottle
ch4 = pyg.headspace_concentration("CH4", 6.85, 1.85, 26.1, 26.2, 1014.7, 35, 140)

# with salinity correction of the Henry's coefficient
co2 = pyg.headspace_concentration(
    "CO2", 1506.7, 403.4, 26.1, 26.2, 1014.7, 35, 140, method="Wiesenburg", salt=24.55
)

# d13C of the dissolved gas (headspace after / before shaking: -14.0 / -8.65 permille)
d13co2 = pyg.headspace_d13c("CO2", 1506.7, 403.4, -14.0, -8.65, 26.1, 26.2, 1014.7, 35, 140)
```

`headspace_d13c` assumes no isotopic fractionation between the water and the headspace after shaking.

### Input types

Inputs can be Python or NumPy scalars, NumPy arrays, lists, pandas Series or [Pint](https://pint.readthedocs.io) quantities. All functions return Pint quantities; use `.magnitude` to get the plain values and `.to("<unit>")` to convert them.

### Units

Units are managed with Pint. Values without units are interpreted with the default units below. Values passed as Pint quantities are converted automatically.

| Key | Default unit |
|---|---|
| `WindSpeed` | m/s |
| `Temperature` | degC |
| `Diss_Gas_Concentration` | mmol m^-3 |
| `Gas_Concentration` | ppm |
| `k600`, `kgas`, `k` | m/d |
| `Atm_Pressure` | hPa |
| `Salinity` | PSU |
| `Diff_Flux` | mmol m^-2 d^-1 |
| `Area` | km^2 |
| `height` | m |
| `Volume` | mL |
| `Isotope_Signature` | permille |

Defaults can be changed for the whole package:
```python
from pyghgaq.functions.units import DEFAULT_INPUT_UNITS

DEFAULT_INPUT_UNITS["WindSpeed"] = "km/h"
DEFAULT_INPUT_UNITS["kgas"] = "cm/h"
```
or for a single call, using the function's argument names as keys:
```python
flux = pyg.atm_diff_flux(0, 1, 2, units={"kgas": "m/s"})
k = pyg.k600(5, "VP2013", area=150, units={"area": "ha"})
```
or by passing Pint quantities:
```python
from pyghgaq.functions.units import Q_

k = pyg.k600(Q_(18, "km/h"))
```

## IMPORTANT
This package is still under development.
