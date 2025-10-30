import yaml


def read_constant(filename: str = "./constant.yml") -> dict:
    """Read constant file.

    Parameters
    ----------
    filename : Path to constant file.

    Return
    -----
    conf_file :
        Constant data.
    """
    with open(filename, "r") as file:
        conf_file = yaml.safe_load(file)
    return conf_file
