"""
This file contains functions for generating theoretical MS1 data for ALife
Polymer Data
"""
from .Constants.GlobalChemicalConstants import *
from .Config_files.Depsipeptide_config import *
from .silico_helpers.insilico_helpers import *


def generate_mass_dictionary(
    monomers,
    max_length,
    ms_level,
    min_length=1,
    sequencing=True,
    universal_rxn=True,
    start_tags=None,
    end_tags=None,
    isobaric_targets=None):
    """
    This function takes a list of monomers and outputs a dictionary of all
    possible sequences or compositions that could arise from input monomers
    within the constraints set, along with associated NEUTRAL monoisotopic
    masses.

    Arguments:
        monomers (list) -- list of monomer one letter codes.
        max_length (int) -- maximum sequence length (in monomer units).

    Keyword Arguments:
        min_length (int) -- minimum sequence length (in monomer units).
                            (default: {1}).
        sequencing (bool) -- specifies whether all possible sequences are to
            be enumerated or just compositions. Set to False if you just need 
            to screen for compositions. (default: {True}).
        universal_rxn (bool) -- specifies whether all monomers are universally
            cross-reactive. If this is set to False, reactivity classes must be 
            read from the polymer config file to ensure that any chemically 
            infeasbile sequences are not included in the final mass dictionary
            (default: {True})
        starting_monomers {list} -- list of monomers. If this list is input,
            only sequences beginning with a monomer in this list will be 
            returned. (default: {None}).
        ending_monomers (list) -- list of monomers. If this is input,
            only sequences ending with a monomer in this list will be returned. 
            (default: {None}.
        isobaric_targets (list) -- list of sequences and / or compositions that
            output sequences must be isobaric to. (default: {None}).

    Returns:
        dict -- dictionary of sequences and associated neutral
            monoisotopic masses. Keys = sequence strings, values = masses 
            (float).
    """

    # check to see if input monomers are set as universally cross-reactive
    # if so, sequence generator used is simpler and (probably) faster
    if universal_rxn:
        sequences = generate_all_sequences(
            monomers=monomers,
            max_length=max_length,
            min_length=min_length,
            sequencing=sequencing,
            start_tags=start_tags,
            end_tags=end_tags,
            isobaric_targets=isobaric_targets
        )
    
    # if input monomers are not universally cross-reactive, more complicated
    # function is required to read polymer config file and generate only
    # those sequences that can be produced from appropriately cross-reactive
    # monomers
    elif not universal_rxn:
        sequences = generate_all_sequences_rxn_classes(monomers)
    
    # work out neutral mass of each sequence and add to mass dictionary in
    # format {seq : mass} where seq = sequence, mass = neutral monoisotopic
    # mass
    massdict = {
        sequence: round(find_sequence_mass(
            sequence=sequence,
            ms_level=1), 
            4)
            for sequence in sequences
    }

    return massdict

def add_adducts_ms1_massdict(
    massdict,
    adducts,
    mode='pos',
    min_z=1,
    max_z=None
):
    """
    Takes an MS1 mass dictionary of sequences and NEUTRAL MASSES and adds
    charged adducts to sequence masses, giving a list of m/z values for
    each sequence.

    Args:
        massdict (dict): dictionary of sequence strings and corresponding
                    neutral masses.
        adducts (list): list of adduct string ids - e.g. ['H', 'Na'].
        mode (str, optional): specifies whether ions with be cationic or
                    anionic. Use 'pos' and 'neg' for positive and negative
                    modes, respectively. Defaults to 'pos'.
        min_z (int, optional): minimum charge of adduct species. Defaults to 1.
        max_z (int, optional): maximum charge of adduct species; if None,
                    max charge is determined by maximum charge of
                    adducts specified in GlobalChemicalConstants.
                    Defaults to None.

    Returns:
        adduct_massdict: dictionary of sequences and corresponding m/z values
                    of charged adduct species.
    """

    # initiate adduct massdict to be populated with sequences and list of m/z
    # values with adducts
    adduct_massdict = {}

    # iterate through sequences, neutral masses and add adducts to
    # adduct_massdict
    for sequence, neutral_mass in massdict.items():
        adduct_massdict[sequence] = add_adducts_sequence_mass(
            mass=neutral_mass,
            adducts=adducts,
            min_z=min_z,
            max_z=max_z,
            mode=mode
        )

    # return ms1 massdict with charged adducts
    return adduct_massdict

def add_loss_products_ms1_massdict(
    massdict,
    max_total_losses=None
):
    """
    This function takes a dictionary of sequences and corresponding NEUTRAL
    MASSES, and incorporates side chain-specific neutral loss products as
    specified in LOSS_PRODUCTS dict in polymer config file.

    Args:
        massdict (dict): dictionary of sequences and corresponding neutral
                        masses.
        max_total_losses (int, optional): maximum number of loss products to
                        incorporate per sequence; if None, all possible loss
                        products for a given sequence are incorporated into
                        output loss_product_dict MS1 mass dictionary.
                        Defaults to None.

    Returns:
        loss_product_dict: dictionary of sequences and neutral mass lists with
                        associated side chain-specific neutral loss products.
    """

    # initiates dictionary to store sequences with associated neutral loss
    # products
    loss_product_dict = {}

    # iterate through sequences and corresponding neutral masses, adding
    # side chain-specific neutral losses to sequence mass lists
    for sequence, neutral_mass in massdict.items():
        loss_product_dict[sequence] = add_sidechain_neutral_loss_products_sequence(
            sequence=sequence,
            sequence_masses=neutral_mass,
            max_total_losses=max_total_losses
        )

    # return dictionary of sequences and associated lists of neutral masses
    # with neutral loss products
    return loss_product_dict

def add_terminal_modification_MS1_sequence(
    sequence,
    sequence_masses,
    terminal_modifications_dict
):
    """
    Takes a sequence, associated masses, modification three letter code and
    associated mass, massdiff, and returns a dict of modified sequence +
    modified sequence masses.

    Args:
        sequence (str): sequence string comprised of monomer one letter codes
        sequence_masses (list): list of monoisotopic NEUTRAL masses for 
            input sequence. 
        terminal_modifications_dict (dict): dict of termini and modifications
    Returns:
        {mod_seq: [mod_seq_masses]}: dictionary of modified sequence string
                        and corresponding modified sequence masses.
    """
    # return nothing if no terminal modifications are specified.
    if not terminal_modifications_dict: 
        return {}
    
    mod_sequences = {}
    
    # iterate through termini in terminal_modifications_dict
    for terminus, terminal_modifications in terminal_modifications_dict.items():
        
        # check whether modifications are specified for current terminus
        if terminal_modifications:

            # retrieve modification information for each modification
            for terminal_modification in terminal_modifications:
                modification_mass=MODIFICATIONS[terminal_modification]["mass"]
                modification_massdiff=MODIFICATIONS[terminal_modification]["mass_diff"]["ms1"]
                modification_termini=MODIFICATIONS[terminal_modification]["termini"]

                # check whether current modification is compatible with 
                # current terminus, as defined MODIFICATIONS in polymer-
                # specific config file
                if float(terminus) not in [float(mod_term) for mod_term in modification_termini]:

                    print(f'you are trying to modify terminus {terminus}')
                    print(f'with {terminal_modification}, but this can only')
                    print(f'modify the following termini: {modification_termini}')
                    raise Exception("please try again")

                # initialise modified sequence dictionary
                mod_sequences = {}

                # modify sequence string
                if type(sequence_masses) != list:
                    sequence_masses = [sequence_masses]
                
                # build terminally modified sequence string
                if float(terminus)==0:
                    mod_seq = f"[{terminal_modification}]{sequence}"

                elif float(terminus)==-1:
                    mod_seq = f"{sequence}[{terminal_modification}]"
                
                # add terminal modification mass to unmodified sequence masses
                modified_sequence_masses = add_modification_sequence_mass_list(
                    mass_list=sequence_masses,
                    modification_mass=modification_mass,
                    mod_mass_diff=modification_massdiff,
                    universal_shift=True
                )
                
                # add modified sequence and the corresponding mass to dictionary
                mod_sequences[mod_seq] = modified_sequence_masses
    
    return mod_sequences
        

def add_terminal_modification_MS1_massdict(
    massdict,
    terminal_modifications_dict
):
    """
    Takes an MS1 sequence mass dictionary and adds a terminal modification to
    every sequence, returning a dictionary of modified sequence string
    and corresponding modified masses.

    Args:
        massdict (dict): dictionary of sequences and corresponding MS1 masses.
        modification (str): modification three letter code from
                MODIFICATIONS_DICT in modificatiosn config file
        modification_mass (float): neutral monoisotopic mass of modifier
        modification_massdiff (float): mass lost upon addition of modification
                to sequence (e.g. H2O for fatty acid acylation of peptides
                as this is a condensation reaction).
        modification_terminus (int): either 0 or -1 to specify start and
                end terminus, respectively.

    Returns:
        modified_massdict (dict): dictionary of modified sequence strings and
                masses.
    """
    #  initialise modified mass dictionary
    modified_massdict = {}

    for sequence, masses in massdict.items():
        modified_massdict.update(
            add_terminal_modification_MS1_sequence(
                sequence=sequence,
                sequence_masses=masses,
                terminal_modifications_dict=terminal_modifications_dict
                )
            )
    
    return modified_massdict

def generate_ms1_mass_dictionary_adducts_losses(
    monomers,
    max_length,
    adducts,
    monomer_modifications_dict,
    universal_monomer_modification,
    universal_terminal_modification,
    terminal_modifications_dict,
    mode='pos',
    min_z=1,
    max_z=None,
    losses=True,
    max_total_losses=None,
    min_length=1,
    chain_terminators=None,
    universal_rxn=True,
    start_tags=None,
    end_tags=None,
    sequencing=True,
    isobaric_targets=None
):
    """ This function generates a loss product dictionary and then
        adds adduct to loss product sequence masses.
    
    Arguments:
        monomers (list) -- list of monomer one letter codes.
        max_length (int) -- maximum sequence length (in monomer units).
        adducts (list): list of adduct string ids - e.g. ['H', 'Na'].
    
    Keyword Arguments:
        mode (str, optional): specifies whether ions with be cationic or
                    anionic. Use 'pos' and 'neg' for positive and negative
                    modes, respectively. Defaults to 'pos'.
        min_z (int, optional): minimum charge of adduct species. Defaults to 1.
        max_z (int, optional): maximum charge of adduct species; if None,
                    max charge is determined by maximum charge of
                    adducts specified in GlobalChemicalConstants.
                    Defaults to None.
        losses (bool) -- Specify whether their is expected to be loss products
                    seen for certain sidechains in the experiments. (default: {True})
        max_total_losses (int) -- maximum number of loss products possible for any given product.
                    (default: {None}).
        loss_product_adducts (dict) -- dictionary of adducts associated with specific loss products.
                    (default: {None}).
        min_length (int) -- minimum length of the product. (default: {1}).
        chain_terminators {[type]} -- [description] (default: {None}).
        universal_rxn (bool) -- [description] (default: {True}).
        start_tags (list) -- list of tags seen at the start of the sequence. (default: {None}).
        end_tags (list) -- list of tags seen at the end of the sequence. (default: {None}).
        sequencing (bool) -- True if you want to attempt to distinguish between different isobarics.
                    If composition is only necessary, set to false. (default: {True}).
        isobaric_targets (list) -- list of specific sequences that you want to or expect to find within the
                    data. (default: {None}).
    
    Returns:
        dict -- dictionary of loss products masses with specified adducts.
    """
    
    # generate list of modified and unmodified monomers 
    monomers = generate_monomers_list_with_sidechain_mods(
        standard_monomers=monomers,
        sidechain_mods=monomer_modifications_dict,
        universal_modification=universal_monomer_modification
    )
    
    # generate neutral mass dictionary of all possible sequences arising from
    # input monomers and constraints set
    MS1_neutral = generate_mass_dictionary(
        monomers=monomers,
        max_length=max_length,
        min_length=min_length,
        ms_level=1,
        sequencing=sequencing,
        universal_rxn=universal_rxn,
        start_tags=start_tags,
        end_tags=end_tags,
        isobaric_targets=isobaric_targets
    )
    
    # check for terminal modificatios; if specified, add these to massdict
    if terminal_modifications_dict.values():

        modified_seq_massdict = add_terminal_modification_MS1_massdict(
            massdict=MS1_neutral,
            terminal_modifications_dict=terminal_modifications_dict
        )

        if universal_terminal_modification:
            if modified_seq_massdict:
                MS1_neutral = modified_seq_massdict
        else:
            MS1_neutral.update(modified_seq_massdict)
    
    # update massdict to include sidechain-specific neutral loss products 
    if losses:
        MS1_neutral = add_loss_products_ms1_massdict(
            massdict=MS1_neutral,
            max_total_losses=max_total_losses
        )
    else:
        MS1_neutral = MS1_neutral

    
    MS1_adduct_dict = add_adducts_ms1_massdict(
        massdict=MS1_neutral,
        adducts=adducts,
        mode=mode,
        min_z=min_z,
        max_z=max_z

    )
    return MS1_adduct_dict

