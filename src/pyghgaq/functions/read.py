import yaml

from importlib import resources
import yaml  # requires pyyaml installed


def read_constant(filename: str = "constant.yml") -> dict:
    """Read constant file.

    Parameters
    ----------
    filename : Path to constant file.

    Return
    -----
    conf_file :
        Constant data.
    """
    with resources.files("pyghgaq").joinpath(filename).open("r") as f:
        constant = yaml.safe_load(f)
    # with open(filename, "r") as file:
    #     conf_file = yaml.safe_load(file)
    return constant
