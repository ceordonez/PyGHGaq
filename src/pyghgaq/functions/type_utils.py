from __future__ import annotations

from typing import TypeAlias, TYPE_CHECKING

from pint import Quantity
from pint.facets.plain import PlainQuantity
from numpy.typing import NDArray

if TYPE_CHECKING:
    pass

Numeric: TypeAlias = float | NDArray | int | Quantity | PlainQuantity
