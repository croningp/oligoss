"""
This file contains functions for generating theoretical MS and MSMS data for
ALife polymer chemistry.
"""
from .MS1_silico import *
from .MS2_silico import *

def __init__():
    """
    Initialises SilicoGenerator and imports constants from polymer config file

    """

def generate_MSMS_sequence_dict(
    monomers,
    max_length,
    fragments,
    adducts,
    mode='pos',
    min_z=1,
    max_z=None,
    losses=True,
    max_neutral_losses=None,
    loss_product_adducts=None,
    min_length=1,
    chain_terminators=None,
    starting_monomers=None,
    ending_monomers=None,
    universal_rxn=True,
    add_signatures=True,
    ms2_signatures=None
):
    """
    This function generates a full MS1-MS2 sequence library

    Args:
        monomers ([type]): [description]
        max_length ([type]): [description]
        fragments ([type]): [description]
        adducts ([type]): [description]
        mode (str, optional): [description]. Defaults to 'pos'.
        min_z (int, optional): [description]. Defaults to 1.
        max_z ([type], optional): [description]. Defaults to None.
        max_neutral_losses (bool, optional): [description]. Defaults to True.
        loss_product_adducts ([type], optional): [description]. Defaults to None.
        min_length (int, optional): [description]. Defaults to 1.
        chain_terminators ([type], optional): [description]. Defaults to None.
        starting_monomers ([type], optional): [description]. Defaults to None.
        ending_monomers ([type], optional): [description]. Defaults to None.
        universal_rxn (bool, optional): [description]. Defaults to True.
        ms2_signatures ([type], optional): [description]. Defaults to None.
    """
    MS1 = generate_ms1_mass_dictionary_adducts_losses(
        monomers,
        max_length,
        adducts,
        mode,
        min_z,
        max_z,
        losses,
        max_neutral_losses,
        loss_product_adducts,
        min_length,
        chain_terminators,
        starting_monomers,
        ending_monomers,
        universal_rxn
    )

    MS2 = generate_ms2_mass_dictionary(
            MS1.keys(),
            fragments,
            adducts,
            mode,
            add_signatures,
            ms2_signatures,
            losses,
            max_neutral_losses,
            loss_product_adducts
    )
    full_MSMS_dict = {}

    for sequence in MS1:
        full_MSMS_dict[sequence] = {
            'MS1': MS1[sequence],
            'MS2': MS2[sequence]
            }
    print(full_MSMS_dict)
