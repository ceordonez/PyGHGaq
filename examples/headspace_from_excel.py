"""Dissolved CH4 and CO2 concentrations and d13C from headspace samples in Excel.

Reads one row per sample from an Excel file, calculates the water
concentration (mmol m^-3) and d13C (permille) of each gas with the
headspace method, and writes the results to a new Excel file.

Requirements (besides pyghgaq): openpyxl
    pip install openpyxl

Usage:
    python headspace_from_excel.py
    python headspace_from_excel.py my_samples.xlsx my_results.xlsx

Edit the settings below to match your columns and bottle volumes.
"""

import sys
from pathlib import Path

import pandas as pd

import pyghgaq as pyg

HERE = Path(__file__).parent

# ---------------------------------------------------------------- settings
INPUT_FILE = HERE / "data" / "headspace_example.xlsx"
OUTPUT_FILE = HERE / "output" / "headspace_example_results.xlsx"

VOL_HS = 35  # headspace volume (mL)
VOL_TOTAL = 140  # total volume of the bottle/syringe (mL)

# Columns common to all gases
TEMP_EQ = "Te_degC"  # equilibration (water) temperature (degC)
TEMP_AIR = "AirT_degC"  # air temperature when the headspace was filled (degC)
PRESSURE = "AirP_hPa"  # atmospheric pressure (hPa)
SALINITY = "sal_psu"  # salinity (PSU); set to None for freshwater

# Columns for each gas:
#   c_hs     headspace concentration after shaking (ppm)
#   c_air    headspace concentration before shaking (ppm)
#   d13c_hs  d13C of the headspace after shaking (permille), or None
#   d13c_air d13C of the headspace before shaking (permille), or None
GASES = {
    "CH4": {
        "c_hs": "HP_12CH4_dry_mean",
        "c_air": "CH4atm_ppm",
        "d13c_hs": "HP_Delta_iCH4_Raw_mean",
        "d13c_air": "d13C_CH4atm_permil",
    },
    "CO2": {
        "c_hs": "12CO2_mean",
        "c_air": "CO2atm_ppm",
        "d13c_hs": "Delta_Raw_iCO2_mean",
        "d13c_air": "d13C_CO2atm_permil",
    },
}
# -------------------------------------------------------------------------


def process(data: pd.DataFrame) -> pd.DataFrame:
    """Add concentration and d13C columns for each gas in GASES."""
    results = data.copy()

    # Henry's coefficient: Weiss 1974 (CO2) and Wiesenburg and Guinasso 1979
    # (CH4) with salinity correction, or Sander 2023 for freshwater.
    if SALINITY is None:
        henry = {"method": "Sanders"}
    else:
        henry = {"method": "Wiesenburg", "salt": data[SALINITY]}

    common = {
        "temp_eq": data[TEMP_EQ],
        "temp_air": data[TEMP_AIR],
        "atmpress": data[PRESSURE],
        "vol_hs": VOL_HS,
        "vol_total": VOL_TOTAL,
        **henry,
    }

    for gas, cols in GASES.items():
        c_hs, c_air = data[cols["c_hs"]], data[cols["c_air"]]

        conc = pyg.headspace_concentration(gas, c_hs, c_air, **common)
        results[f"{gas}_mmolm3"] = conc.to("mmol m^-3").magnitude

        if cols.get("d13c_hs") and cols.get("d13c_air"):
            d13c = pyg.headspace_d13c(
                gas, c_hs, c_air, data[cols["d13c_hs"]], data[cols["d13c_air"]],
                **common,
            )
            results[f"d13C_{gas}_permil"] = d13c.to("permille").magnitude

    return results


def main() -> None:
    input_file = Path(sys.argv[1]) if len(sys.argv) > 1 else INPUT_FILE
    output_file = Path(sys.argv[2]) if len(sys.argv) > 2 else OUTPUT_FILE

    data = pd.read_excel(input_file)
    results = process(data)

    output_file.parent.mkdir(parents=True, exist_ok=True)
    results.to_excel(output_file, index=False, sheet_name="Results")
    print(f"Results written to {output_file}")
    new_cols = [c for c in results.columns if c not in data.columns]
    print(results[[data.columns[0], *new_cols]].to_string(index=False))


if __name__ == "__main__":
    main()
