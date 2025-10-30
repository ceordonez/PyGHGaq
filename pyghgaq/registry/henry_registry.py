from functools import wraps
from typing import Any, Callable

import numpy as np

type exporterfunc = Callable[..., np.ndarray | float]
type exporterregis = dict[str, dict[str, Callable]]
exporters: exporterregis = {}


def register_hcp(varname: str, salt: str):
    def decorator(func: Callable[..., Any]):
        @wraps(func)
        def wrapper(*args: Any, **kwards: Any) -> Any:
            return func(*args, **kwards)

        exporters.setdefault(varname, {})[salt] = wrapper
        return wrapper

    return decorator
