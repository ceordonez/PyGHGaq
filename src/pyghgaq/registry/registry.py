import inspect
from copy import deepcopy
from functools import wraps
from typing import Any, Callable

import numpy as np

from pyghgaq.functions.type_utils import Numeric
from pyghgaq.functions.units import Q_, to_si

type exporterfunc = Callable[..., np.ndarray | float]
type exporterhcp = dict[str, dict[str, Callable]]
type exporterk600 = dict[str, exporterfunc]
exportershcp: exporterhcp = {}
exportersk600: exporterk600 = {}
exporterssh: exporterk600 = {}


def register_hcp(varname: str, salt: str):
    def decorator(func: Callable[..., Any]):
        @wraps(func)
        def wrapper(*args, **kwargs):
            return func(*args, **kwargs)

        exportershcp.setdefault(varname, {})[salt] = wrapper
        return func

    return decorator


def register_k600(method: str):
    def decorator(func: Callable[..., Any]):
        @wraps(func)
        def wrapper(*args: Any, **kwargs: Any) -> Any:
            return func(*args, **kwargs)

        exportersk600[method] = wrapper
        return wrapper

    return decorator


def register_schmidt(method: str):
    def decorator(func: Callable[..., Any]):
        @wraps(func)
        def wrapper(*args: Any, **kwargs: Any) -> Any:
            return func(*args, **kwargs)

        exporterssh[method] = wrapper
        return wrapper

    return decorator


def enforce_units(**expected_units):
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            from inspect import signature

            sig = signature(func)
            bound = sig.bind(*args, **kwargs)
            bound.apply_defaults()

            for name, exp_unit in expected_units.items():
                if name in bound.arguments and isinstance(
                    bound.arguments[name], Numeric
                ):
                    bound.arguments[name] = Q_(
                        bound.arguments[name], bound.arguments["units"][name]
                    )
                    bound.arguments[name] = to_si(bound.arguments[name], exp_unit)
            return func(*bound.args, **bound.kwargs)

        return wrapper

    return decorator


def warn_default_units(func):
    @wraps(func)
    def wrapper(*args, units=None, **kwargs):
        # Retrieve the default units from the function signature
        sig = inspect.signature(func)
        param = sig.parameters.get("units")
        if param is None or param.default is None:
            # No default units defined in the function
            default_units = {}
        else:
            default_units = deepcopy(param.default)  # copy to avoid mutation

        if units is None:
            units = {}

        # Warn for each default unit not explicitly provided by the user
        for key, default in default_units.items():
            if key not in units:
                print(f"Assuming {key} = '{default}' (default)")

        # Merge user-provided units with defaults
        merged_units = {**default_units, **units}

        # Call the original function with the merged units
        return func(*args, units=merged_units, **kwargs)

    return wrapper
