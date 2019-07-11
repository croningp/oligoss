"""
This file contains functions for generating theoretical MS and MSMS data for 
ALife polymer chemistry. 
"""
from GlobalChemicalConstants import *
from Config_files.Depsipeptide_config import *


def __init__():
    """
    Initialises SilicoGenerator and imports constants from polymer config file

    """
    #ADD SOMETHING TO IMPORT CONFIG FILE HERE - IMPORTING DEPSIPEPTIDE 
    #CONFIG FILE DIRECTLY IS A TEMPORARY SOLUTION FOR INITIAL TESTING OF 
    #FUNCTIONS
    
    pass

def enumerate_all_sequences(
    monomers,
    max_length,
    min_length=1, 
    sequencing=True,
    chain_terminators=None,
    starting_monomers=None,
    ending_monomers=None):
    """
    [This function takes a list of input monomers and outputs all possible
    sequences or compositions that could arise within the constraints set, 
    which are described below]
    
    Arguments:
        monomers {list} -- [list of monomer one letter codes]
        max_length {int} -- [maximum sequence length (in monomer units)]
    Keyword Arguments:
        min_length {int} -- [minimum sequence length (in monomer units)]
                            (default: {1})
        sequencing {bool} -- [specifies whether all possible sequences are to 
                            be enumerated or just compositions. Set to False if 
                            you just need to screen for compositions] 
                            (default: {True})
        chain_terminators {list} -- [list of monomers which prevent further
                            elongation] (default: {None})
        starting_monomers {list} -- [list of monomers. If this list is input,
                            only sequences beginning with a monomer in this 
                            list will be returned] (default: {None})
        ending_monomers {list} -- [list of monomers. If this is input, 
                            only sequences ending with a monomer in this list
                            will be returned] (default: {None})
    """
    sequences = monomers
    
    for i in range(0, max_length):
        new_sequences = []
        for monomer in monomers: 
            new_sequences.extend([monomer+sequence for sequence in sequences])
        sequences.extend(new_sequences)
        sequences = list(set(sequences))

    if chain_terminators: 
        sequences = [
            seq for seq in sequences 
            if (seq[0] in chain_terminators
            or seq[-1] in chain_terminators)
        ]
    if starting_monomers:
        sequences = [
            seq for seq in sequences 
            if seq[0] in starting_monomers
        ]
    if ending_monomers: 
        sequences = [
            seq for seq in sequences 
            if seq[-1] in ending_monomers
        ]

    if not sequencing:   
        sequences = list(set(["".join(sorted(seq)) for seq in sequences]))

    sequences = [seq for seq in sequences if len(seq) >= min_length]
    
    return sequences 

def find_sequence_mass(sequence):
    """
    [Takes a sequence string and returns neutral monoisotopic mass for sequence]
    
    Arguments:
        sequence {str} -- [sequence string made up of constituent monomer one 
                            letter codes]
    
    Returns:
        sequence_mass {float} -- [neutral monoisotopic mass of input sequence]
    """
    sequence_mass = sum([MONOMERS[c][0] for i, c in enumerate(sequence)])
    sequence_mass -= (len(sequence)-1)*MASS_DIFF
    return sequence_mass

def enumerate_all_sequences_rxn_classes(monomers):
    """
    [This function will be written to build sequence lists in cases where 
    monomers are not universally cross-reactive - i.e. more complicated 
    operations are required]
    
    Arguments:
        monomers {[type]} -- [description]
    """
    sequences = []
    return sequences

def generate_mass_dictionary(
    monomers,
    max_length,
    min_length=1, 
    sequencing=True,
    universal_rxn=True,
    chain_terminators=None,
    starting_monomers=None,
    ending_monomers=None): 
    """
    [This function takes a list of monomers and outputs a dictionary of all 
    possible sequences or compositions that could arise from input monomers
    within the constraints set, along with associated NEUTRAL monoisotopic 
    masses]
    
    Arguments:
        monomers {list} -- [list of monomer one letter codes]
        max_length {int} -- [maximum sequence length (in monomer units)]
    
    Keyword Arguments:
        min_length {int} -- [minimum sequence length (in monomer units)]
                            (default: {1})
        sequencing {bool} -- [specifies whether all possible sequences are to 
                            be enumerated or just compositions. Set to False if 
                            you just need to screen for compositions] 
                            (default: {True})
        universal_rxn {bool} -- [specifies whether all monomers are universally
                            cross-reactive. If this is set to False, reactivity
                            classes must be read from the polymer config file
                            to ensure that any chemically infeasbile sequences
                            are not included in the final mass dictionary]
                            (default: {True})
        chain_terminators {list} -- [list of monomers which prevent further
                            elongation] (default: {None})
        starting_monomers {list} -- [list of monomers. If this list is input,
                            only sequences beginning with a monomer in this 
                            list will be returned] (default: {None})
        ending_monomers {list} -- [list of monomers. If this is input, 
                            only sequences ending with a monomer in this list
                            will be returned] (default: {None})
    
    Returns:
        [massdict] -- [dictionary of sequences and associated neutral 
                    monoisotopic masses. Keys = sequence strings, values = 
                    masses (float)]
    """

    if not universal_rxn: 
        sequences = enumerate_all_sequences(
            monomers,
            max_length,
            min_length,
            sequencing,
            chain_terminators,
            starting_monomers,
            ending_monomers
            )
    elif universal_rxn: 
        sequences = enumerate_all_sequences_rxn_classes(monomers)
    
    massdict = dict(
        (sequence, find_sequence_mass(sequence)) 
        for sequence in sequences
        )
    
    return massdict 


"""
[INSERT LOSS PRODUCT FUNCTION(S) HERE]
"""
"""
[FUNCTIONS TO PRODUCE READING FRAME SHIFTS FOR SEQUENCES] 

"""

def build_fragment_series_single_sequence(
    sequence, 
    fragment_series, 
    mode='pos'):
    """
    [This function takes a sequence and builds an MS2 fragment series,
    returning a dictionary of fragment ids and corresponding m/z values]
    
    Arguments:
        sequence {str} -- [sequence string of constituent monomer one letter
                            codes]
        fragment_series {str} -- [fragment string used to denote fragment 
                                series]

    Keyword Arguments:
        mode {str} -- [specifies whether fragmented species are cationic or
                    anionic, use 'pos' for cations and 'neg' for anions] 
                    (default: {'pos'})
    
    Returns:
        [fragment_dict] -- [dictionary of fragments and corresponding m/z 
                            values]
    """
    terminus = FRAG_SERIES[fragment_series]['terminus']
    mass_diff = FRAG_SERIES[fragment_series]['mass_diff']['mode']
    start = FRAG_SERIES[fragment_series]['start']
    end = FRAG_SERIES[fragment_series]['end']
    
    if terminus == -1:
        sequence = sequence[::-1]
    
    #build a fragment dictionary for sequence
    fragment_dict = {
        f'{fragment_series}{i+1}':
        find_sequence_mass(sequence[0:i]) + mass_diff
        for i in range(0 + start, len(sequence) - end)
    } 
    
    return fragment_dict 

def generate_monomer_signature_ions_sequence(
    sequence, 
    signature,
):
    """
    [This function generates a dictionary of monomer signature ions for 
    a sequence]
    
    Arguments:
        sequence {str} -- [sequence string with constituent monomer one letter
                            codes]
        signature {str} -- [signature fragment code - e.g. 'Im' for amino acid
                            immonium fragments]
    
    Returns:
        [MS2_signature_dict] -- [dictionary of monomer signature fragments
                            and associated m/z values]
    """
    monomers = list(set([c for i, c in enumerate(sequence)]))

    signature_ions = MS2_SIGNATURE_IONS[signature]

  
    MS2_signature_dict = {
        f'{signature}{monomer}': signature_ions[monomer]
        for monomer in monomers
        if monomer in signature_ions
    }
    return MS2_signature_dict