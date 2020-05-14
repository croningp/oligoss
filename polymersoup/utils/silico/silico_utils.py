"""
This file contains utilities for the silico module
"""
from ..global_chemical_constants import CATIONS, ANIONS

#  dict of modification markers for denoting sidechain, terminal and crosslink-
#  modified monomers in sequence strings
SILICO_MOD_MARKERS = {
    '(': ')',  # marks beginning and end of sidechain modification
    'terminal': ['[', ']']  # marks terminal modification
}


def retrieve_adduct_info(mode, adducts):
    """
    Takes a list of adduct strings, Parameters object and returns dict of
    adduct strings and adduct info: mass, charge range.

    Args:
        mode (str): either "pos" or "neg" for positive and negative mode,
            respectively
        adducts (List[str]): list of adduct strings

    Raises:
        Exception: raised if params.mode not in ["pos", "neg"]

    Returns:
        Dict[str, List[float]]: dict of adduct strings and associated info.
            info = [mass (float), minimum_charge (int), maximum_charge (int)]
    """
    if mode == "pos":
        ions = CATIONS
    elif mode == "neg":
        ions = ANIONS
    else:
        raise Exception(
            f'mode must be either "pos" or "neg", {mode} is not an\n'
            'option!'
        )

    return {
        adduct: ions[adduct]
        for adduct in adducts
    }
