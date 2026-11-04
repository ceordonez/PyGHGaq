from __future__ import annotations

from typing import TYPE_CHECKING, TypeAlias

import numpy as np
from pint import Quantity
from pint.facets.plain import PlainQuantity

if TYPE_CHECKING:
    pass

Numeric: TypeAlias = float | np.ndarray | int | Quantity | PlainQuantity
