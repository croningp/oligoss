"""
This file contains essential functions which are essential to all
in silico operations
"""
from ..Constants.GlobalChemicalConstants import *
from ..Config_files.Depsipeptide_config import *
import copy

def find_sequence_mass(
    sequence,
    four_point_float=True
):
    """
    [Takes a sequence string and returns neutral monoisotopic mass for sequence]

    Arguments:
        sequence {str} -- [sequence string made up of constituent monomer one
                            letter codes]
        decimal_points {int} -- specifies number of decimal points to be worked
                            out for sequence mass
    Returns:
        sequence_mass {float} -- [neutral monoisotopic mass of input sequence]
    """
    sequence_mass = sum([MONOMERS[c][0] for i, c in enumerate(sequence)])
    sequence_mass -= (len(sequence)-1)*MASS_DIFF
    if four_point_float:
        sequence_mass = float(f"{sequence_mass:.4f}")
    return sequence_mass

def add_adducts_sequence_mass(
    neutral_mass,
    adducts,
    min_z=1,
    max_z=None,
    mode='pos'
):
    """
    [This functions adds charged adducts to a neutral sequence mass]

    Arguments:
        neutral_mass {float} -- [neutral monoisotopic mass of sequence]
        adducts {list} -- [list of adduct strings. All adducts must be found
                        in either ANIONS or CATIONS dicts in
                        GlobalChemicalConstants.py]
    Keyword Arguments:
        min_z {int} -- [minimum ABSOLUTE charge of adduct] (default: {1})
        max_z {int} -- [maximum ABSOLUTE charge of adduct. If this is set to
                        None, maximum charge is assumed to be maximum oxidation
                        state of adduct ions] (default: {None})
        mode {str} -- [either 'pos' or 'neg' for positive and negative mode
                        mass spectrometry, respectively] (default: {'pos'})
    Returns:
        [charged_sequence_masses] -- [list of m/z values for charged sequences
                        with adducts]
    """
    #list anions and cations within adducts
    anions = [adduct for adduct in adducts if adduct in ANIONS]
    cations = [adduct for adduct in adducts if adduct in CATIONS]

    if type(neutral_mass) != list:
        neutral_mass = [neutral_mass]

    #check whether counterions need to be considered for adducts which have the
    #opposite charge from mode. If so, return charged adducts from adduct_complex
    #function
    if (mode == 'pos' and len(anions) > 0) or (mode == 'neg' and len(cations) > 0):
        charged_sequence_masses = add_adduct_complexes_sequence_mass(
            sequence,
            neutral_mass,
            adducts,
            min_z,
            max_z)
        return charged_sequence_masses

    #retrieve adduct masses and charges from GlobalChemicalConstants
    if mode == 'pos':
        ions = CATIONS
    elif mode == 'neg':
        ions = ANIONS

    #initiate list to add charged adduct m/z values
    charged_sequence_masses = []

    #iterate through adducts and return m/z values of charged species
    #include multiply charged species that fit within constrains of min_z, max_z
    #and ion oxidation states given in GlobalChemicalConstants
    for adduct in adducts:
        adduct_mass = ions[adduct][0]
        min_charge = ions[adduct][1]
        max_charge = ions[adduct][2]

        # get maximum and minimum charge - either specified or use default
        # minimum and maximum charge states of adducts from
        # GlobalChemicalConstants
        if not max_z:
            max_z = max_charge
        min_z, max_z = max([min_z, min_charge]), min([max_z, max_charge])

        for i in range(min_z, max_z+1):
            charged_sequence_masses.extend(
                [(mass+adduct_mass)/i for mass in neutral_mass])

    #finally, return list of charged sequence m/z values
    return sorted(list(set(
                [float(f"{mass:0.4f}")
                for mass in charged_sequence_masses])))

def add_adduct_complexes_sequence_mass(
    sequence,
    neutral_mass,
    adducts,
    min_z=1,
    max_z=None,
    mode='pos'
    ):
    """
    [This function will add complex adducts to sequence masses - i.e. multiple
    metal centres, ions with counterions and associated solvent species etc]

    Arguments:
        sequence {[type]} -- [description]
        neutral_mass {[type]} -- [description]
        adducts {[type]} -- [description]

    Keyword Arguments:
        min_z {int} -- [description] (default: {1})
        max_z {[type]} -- [description] (default: {None})
        mode {str} -- [description] (default: {'pos'})

    Returns:
        [type] -- [description]
    """
    raise Exception('this function is not complete - do NOT USE')
    return neutral_mass

def generate_all_sequences(
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
    Returns:
        sequences {list} -- [list of possible sequences that could arise from
                            input monomers within the constraints set]
    """
    sequences = copy.deepcopy(monomers)
    for i in range(0,max_length-1):
        for monomer in monomers:
            sequences.extend(
                [seq+monomer for seq in sequences
                if len(seq) == i +1])
    if not sequencing:
        sequences = ["".join(sorted(seq)) for seq in sequences]
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

def generate_all_sequences_rxn_classes(monomers):
    """
    [This function will be written to build sequence lists in cases where
    monomers are not universally cross-reactive - i.e. more complicated
    operations are required]

    Arguments:
        monomers {[type]} -- [description]
    """
    sequences = []
    return sequences

def add_sidechain_neutral_loss_products_sequence(
    sequence,
    sequence_masses,
    max_total_losses=None
):
    """
    [Takes a sequence, mass and subtracts monomer-specific side chain neutral
    loss products as specified in LOSS_PRODUCTS dictionary found in polymer
    config file]

    Arguments:
        sequence {str} -- [sequence string comprised of monomer one letter
                            codes]
        sequence_masses {float or list} -- [float for a single mass, or list
                                    of floats for multiple associated masses,
                                    e.g. for various adducts]

    Keyword Arguments:
        max_total_losses {int} -- [specifies the maximum number of loss
                                products to be subtracted from a given sequence
                                if None, all possible loss products will be
                                incorporated into final mass list]
                                (default: {None})

    Returns:
        [final_masses] {list} -- [list of m/z values for sequence +/- side
                                chain neutral losses]
    """
    # make sequence_mass a list if not already
    if type(sequence_masses) != list:
        sequence_masses = [sequence_masses]

    # get dictionary of loss-product prone monomers and associated loss products
    monomers = list(set([c for i, c in enumerate(sequence)]))
    loss_monomers = {
        monomer: LOSS_PRODUCTS[monomer]
        for monomer in monomers
        if monomer in LOSS_PRODUCTS
    }

    # initiate list of final sequence masses, including loss products AND full
    # sequence mass(es)
    final_masses = sequence_masses

    # keep track of number of loss products subtracted from a sequence -
    # required if there is an upper cap on number of loss products specified by
    # max_total_losses argument
    total_losses = 0

    # iterate through loss-product prone monomers and subtract associated loss
    # products
    for monomer in loss_monomers:

        # count occurence of each loss-prone monomer in the sequence
        occurence = sequence.count(monomer)

        # initiate list of side chain loss products (sequence_mass-loss) for
        # monomer
        monomer_losses = []

        # iterate through associated side chain losses for monomer and subtract
        # from sequence masses

        for loss in loss_monomers[monomer]:
            for i in range(1, occurence+1):
                neutral_loss = loss*i
                monomer_losses.extend(
                    [
                        mass-neutral_loss
                        for mass in final_masses
                    ]
            )
                # keep count of total losses subtracted. if this reaches the
                # limit specified by max_total_losses (if specified), stop
                # adding loss products and return final masses list as is
                total_losses += 1
                if max_total_losses and total_losses >= max_total_losses:
                    final_masses.extend(monomer_losses)

                    # make masses 4 point floats
                    final_masses = [
                        float(f"{mass:.4f}")
                        for mass in final_masses
                    ]

                    # return sorted list of loss products, removing duplicates
                    return sorted(list(set(final_masses)))

        # add monomer-specific loss products to final sequence masses,
        # remove duplicate final masses
        # make masses four point floats
        final_masses.extend(monomer_losses)
        final_masses = sorted(
            list(
                set(
                    [float(f"{mass:0.4f}") for mass in final_masses])))

    return final_masses

def generate_reading_frames_sequence(sequence):
    """
    Takes a linear sequence and outputs a list of reading frame shifts for
    later use in generating proposed fragments for cyclic sequences. A reading
    frame shift is carried out by taking the final monomer in the sequence
    and making it the first monomer, shifting the remaining monomers up 1 index
    in the polymer chain. Example: one reading frame shift of sequence 'AAV'
    produces shifted sequence 'VAA'

    Args:
        sequence (str): sequence string consisting of monomer one letter codes

    Returns:
        reading_frames: list of unique reading frames
    """
    # initiate reading frames list with input sequence
    reading_frames = [sequence]

    # get list of unique monomers within sequence; if only kind of one monomer
    # is present, there can only be one reading frame for the sequence.
    # Therefore, sequence is returned
    unique_monomers = list(set([c for i, c in enumerate(sequence)]))
    if len(unique_monomers) == 1:
        return reading_frames

    # one by one, generate reading frames and check if they are unique (i.e.
    # not yet in reading_frames list); if so, add to list of reading frames
    for i in range(0, len(sequence)+1):
        # take last monomer in sequence (sequence[-1]) and make it the first
        # monomer
        sequence = sequence[-1] + sequence[0:len(sequence)-1]
        if sequence not in reading_frames:
            reading_frames.append(sequence)

    # return list of unique reading frames
    return reading_frames

def generate_dict_isobaric_sequences(sequences):
    """
    Takes a list of sequences and groups isobaric sequences in a dictionary.
    Output = dictionary of isobaric groups (key = sorted sequence, value = list
    of sequences isobaric to sorted sequence)

    Args:
        sequences (list): list of sequence strings
    """
    # get list of sequence compositions
    sorted_sequences = ["".join(sorted(sequence)) for sequence in sequences]

    # generate dictionary of sequence compositions (key = sorted sequences)
    # and lists of sequences matching those compositions (value = sequence list)
    isobaric_dict = {
        sorted_seq: [
            seq for seq in sequences
            if "".join(sorted(seq)) == sorted_seq
        ]
            for sorted_seq in sorted_sequences
        }
    return isobaric_dict
