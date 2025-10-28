# PyGHGaq 

This package contains different methods and properties to estimate GHG (CH<sub>4</sub> and CO<sub>2</sub>) emissions and concentrations for aquatic ecosystems

## Installation
```sh
pip install pyghgaq
```

## Usage

```python
import pyghgaq
import numpy as np

cw = np.random.random(20)
csat = np.random.random(20)
k = np.random.random(20)
flux = pyghgaq.atm_diff_flux(cw, csat, k)
```

## IMPORTANT
This package is still under delvelopment
