import importlib

import numpy as np


def schmidt_number(varname: str, temp_c: np.ndarray | float):
    from pyghgaq.registry.registry import exporterssh

    importlib.import_module(f"pyghgaq.gases.{varname}")
    exporter = exporterssh.get(varname)
    if exporter is None:
        raise ValueError(
            f"Schmit number calcultion for '{varname}' gas has not been implemented"
        )
    return exporter(temp_c)
