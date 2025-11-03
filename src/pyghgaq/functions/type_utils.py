from typing import TypeAlias

import numpy as np
from pint import Quantity

Numeric: TypeAlias = float | np.ndarray | Quantity | int
