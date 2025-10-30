from functools import wraps
from typing import Any, Callable

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
