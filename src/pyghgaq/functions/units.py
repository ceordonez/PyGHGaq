from pint import UnitRegistry

ureg = UnitRegistry(autoconvert_offset_to_baseunit = True)
Q_ = ureg.Quantity

## HERE YOU CAN DEFINE MORE UNITS

ureg.formatter.default_format = '~P'

def to_si(value, expected_unit):
    """
    Ensure `value` has units of `expected_unit`, then convert to SI.
    
    Parameters
    ----------
    value : float | pint.Quantity
        Value with or without units.
    expected_unit : str | pint.Unit
        Expected physical unit (e.g., 'm/s', 'K', 'Pa').
    
    Returns
    -------
    pint.Quantity
        Quantity converted to SI base units.
    """
    expected_unit = ureg.Unit(expected_unit)

    # If value is unitless, assign expected units
    if not isinstance(value, ureg.Quantity):
        q = Q_(value, expected_unit)
    else:
        q = value.to(expected_unit)
    
    __import__('pdb').set_trace()
    # Convert to SI base units
    return q.to_base_units()
