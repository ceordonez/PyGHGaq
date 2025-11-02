from functools import wraps
from typing import Any, Callable
from pyghgaq.functions.units import to_si

import numpy as np

type exporterfunc = Callable[..., np.ndarray | float]
type exporterhcp = dict[str, dict[str, Callable]]
type exporterk600 = dict[str, exporterfunc]
exportershcp: exporterhcp = {}
exportersk600: exporterk600 = {}
exporterssh: exporterk600 = {}


def register_hcp(varname: str, salt: str):
    def decorator(func: Callable[..., Any]):
        @wraps(func)
        def wrapper(*args: Any, **kwards: Any) -> Any:
            return func(*args, **kwards)

        exportershcp.setdefault(varname, {})[salt] = wrapper
        return wrapper

    return decorator


def register_k600(method: str):
    def decorator(func: Callable[..., Any]):
        @wraps(func)
        def wrapper(*args: Any, **kwards: Any) -> Any:
            return func(*args, **kwards)

        exportersk600[method] = wrapper
        return wrapper

    return decorator


def register_schmidt(method: str):
    def decorator(func: Callable[..., Any]):
        @wraps(func)
        def wrapper(*args: Any, **kwards: Any) -> Any:
            return func(*args, **kwards)

        exporterssh[method] = wrapper
        return wrapper

    return decorator

def enforce_units(**expected_units):
    def decorator(func):
        def wrapper(*args, **kwargs):
            from inspect import signature
            sig = signature(func)
            bound = sig.bind(*args, **kwargs)
            bound.apply_defaults()

            for name, exp_unit in expected_units.items():
                if name in bound.arguments:
                    bound.arguments[name] = to_si(bound.arguments[name], exp_unit)
            return func(*bound.args, **bound.kwargs)
        return wrapper
    return decorator
