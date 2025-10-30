import importlib

import numpy as np

# import inspect


# def _update_henry_coefficient_doc(henry_coefficient):
#     """Generate a dynamic docstring listing all registered exporters."""
#     ## NOT IMPLEMENTED YET
#
#     lines = ["Available exporters and their expected parameters:\n"]
#
#     from registry.registry import exportershcp
#
#     for varname, methods in exportershcp.items():
#         for method, func in methods.items():
#             sig = inspect.signature(func)
#             lines.append(f"  {varname} / {method}{sig}")
#     henry_coefficient.__doc__ = "\n".join(lines)
#


def schmidt_number(varname: str, temp_c: np.ndarray | float):
    from pyghgaq.registry.registry import exporterssh

    importlib.import_module(f"pyghgaq.gases.{varname}")
    exporter = exporterssh.get(varname)
    if exporter is None:
        raise ValueError(
            f"Schmit number calcultion for '{varname}' gas has not been implemented"
        )
    return exporter(temp_c)
