"""
This file contains constants and functions relevant to handling parameters
from polymer class-specific config files

"""
import os
import json

HERE = os.path.dirname(__file__)

POLYMER_CONFIGS = os.path.join(
    os.path.dirname(HERE), "polymer_configs")

POLYMER_ALIASES = {
    f.replace(".json", ""): os.path.join(POLYMER_CONFIGS, f)
    for f in os.listdir(POLYMER_CONFIGS)
}

def load_polymer_info(polymer_class):
    """
    Returns polymer_info for polymer_class

    Args:
        polymer_class (str): alias for polymer backbone (i.e. name of file
            in polymer_configs folder)

    Raises:
        Exception: raised if polymer_class not a valid alias

    Returns:
        dict: polymer_info dict
    """
    #  if polymer_class is not a config not alias, return immediately
    if type(polymer_class) == dict:
        return polymer_class

    #  check whether polymer_class str is valid alias for polymer backbone props
    polymer_class = polymer_class.replace(".json", "")
    if polymer_class not in POLYMER_ALIASES:
        raise Exception(
            f'{polymer_class} not a valid polymer alias. Please choose from\n'
            f'the following: {POLYMER_ALIASES.keys()}')

    #  retrieve polymer_info dict and return
    with open(POLYMER_ALIASES[polymer_class]) as f:
        return json.load(f)
