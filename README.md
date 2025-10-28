# PyGHGaq 

This package contains different methods and properties to estimate GHG ($\text{CH}_4$ and $\text{CO}_2$) emissions and concentrations for aquatic ecosystems

## Installation
pip install pyghgaq

## Usage

```python
import pyghgaq

cw = np.random(20)
csat = np.random(20)
k = np.random(20)
flux = pyghgaq.atm_diff_flux(cw, csat, k)
```

## IMPORTANT
THIS PACKAGE IS STILL IN DEVELOPMENT
