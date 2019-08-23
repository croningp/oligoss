"""
This file contains essential functions which are essential to all
in silico operations
"""
from ..Constants.GlobalChemicalConstants import *
from ..Config_files.Depsipeptide_config import *
import copy
import time
import itertools

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

    # calculate sequence mass by summing the mass of its constituent monomers,
    # subtracting the monomer addition MASS_DIFF for every monomer addition
    sequence_mass = sum([MONOMERS[c][0] for i, c in enumerate(sequence)])
    sequence_mass -= (len(sequence)-1)*MASS_DIFF

    # check whether sequence mass is to be trimmed to four decimal places
    if four_point_float:
        sequence_mass = round(sequence_mass, 4)
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
        neutral_mass {float} -- neutral monoisotopic mass of sequence
        adducts {list} -- list of adduct strings. All adducts must be found
                        in either ANIONS or CATIONS dicts in
                        GlobalChemicalConstants.py
    Keyword Arguments:
        min_z {int} -- minimum ABSOLUTE charge of adduct (default: {1})
        max_z {int} -- maximum ABSOLUTE charge of adduct. If this is set to
                        None, maximum charge is assumed to be maximum oxidation
                        state of adduct ions (default: {None})
        mode {str} -- either 'pos' or 'neg' for positive and negative mode
                        mass spectrometry, respectively (default: {'pos'})
    Returns:
        [charged_sequence_masses] -- list of m/z values for charged sequences
                        with adducts
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
        # add more comments
        min_z, max_z = max([min_z, min_charge]), min([max_z, max_charge])

        for i in range(min_z, max_z+1):
            charged_sequence_masses.extend(
                [(mass+adduct_mass)/i for mass in neutral_mass])

    #finally, return list of charged sequence m/z values
    return sorted(list(set(
                [float(f"{mass:0.4f}")
                for mass in charged_sequence_masses])))

def find_functional_groups_monomer(
    monomer
):
    """ This function returns the functional groups present in a monomer

    Arguments:
        monomer {string} -- One letter monomer code that has associated neutral
        monoisotopic mass, list of functional groups with their associated
        frequency per monomer.

    Returns:
        func_groups -- List of functional groups present in the monomer.
    """

    monomer_info = MONOMERS[monomer]

    func_groups = [lambda x: x[0] for x in monomer_info[1]]

    return func_groups

def add_adduct_complexes_sequence_mass(
    sequence,
    neutral_mass,
    adducts,
    min_z=1,
    max_z=None,
    mode='pos',
    max_total_ions=None,
    min_total_ions=None,
    max_mode_matching_ions=None,
    min_mode_matching_ions=None,
    max_opposite_charged_ions=None,
    min_opposite_charged_ions=None,
    max_mz_overall_adduct=0,
    min_mz_overall_adduct=None
    ):
    """
    This function will add complex adducts to sequence masses - i.e. multiple
    metal centres, ions with counterions and associated solvent species etc

    Args:
        sequence (str): sequence string comprised of monomer one letter codes
        neutral_mass (float): neutral monoisotopic mass of sequence
        adducts (list): list of adduct strings
        min_z (int, optional): minimum absolute charge of OVERALL ION COMPLEX.
                    Defaults to 1.
        max_z ([type], optional): maximum absolute charge of OVERALL ION COMPLEX.
                    Defaults to None.
        mode (str, optional): specifies whether OVERALL ION COMPLEX is positive
                    or negative; either 'pos' or 'neg' for positive and
                    negative mode, respectively. Defaults to 'pos'.
        max_total_ions (int, optional): maximum number of total ions of
                    ANY CHARGE STATE IN OVERALL ION COMPLEX. Defaults to None.
        min_total_ions (int, optional): minimum number of total ions of
                    ANY CHARGE STATE IN OVERALL ION COMPLEX. Defaults to None.
        max_mode_matching_ions (int, optional): maximum number of ions whose
                    CHARGE STATE MATCHES THE OVERALL CHARGE STATE OF THE
                    ADDUCT. Defaults to None.
        min_mode_matching_ions (int, optional): minimum number of ions whose
                    CHARGE STATE MATCHES THE OVERALL CHARGE STATE OF THE
                    ADDUCT. Defaults to None.
        max_opposite_charged_ions (int, optional): maximum number of ions
                    whose CHARGE STATE IS OPPOSITE TO THE OVERALL CHARGE STATE
                    OF THE ADDUCT. Defaults to None.
        min_opposite_charged_ions ([type], optional): minimum number of ions
                    whose CHARGE STATE IS OPPOSITE TO THE OVERALL CHARGE STATE
                    OF THE ADDUCT. Defaults to None.

    Raises:
        Exception: [description]

    Returns:
        [type]: [description]
    """

    # initialise list of final m/z values of adduct complexes
    masses = []

    if not max_mz_overall_adduct:
        max_mz_overall_adduct = sorted(masses)[-1]

    # filters out m/z values that do not fall within the specific minimum and
    # maximum m/z range
    masses = filter(lambda mass: mass >= min_mz_overall_adduct and mass <=
                max_mz_overall_adduct, masses)

    return masses

def generate_all_sequences(
    monomers,
    max_length,
    min_length=1,
    sequencing=True,
    chain_terminators=None,
    start_tags=None,
    end_tags=None,
    isobaric_targets=None):
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
        start_tags {list} -- list of monomers. If this list is input, sequences
                            are tagged at terminus 0 by each of the monomers
                            in start_tags - one tag per sequence, and only
                            tagged sequences are returned (default: {None})
        end_tags {list} -- list of monomers. If this list is input, sequences
                            are tagged at terminus -1 by each of the monomers
                            in end_tags - one tag per sequence, and only
                            tagged sequences are returned (default: {None})
        isobaric_targets {list} -- list of sequences and / or compositions that
                            output sequences must be isobaric to. (default: 
                            {None})
    Returns:
        sequences {list} -- [list of possible sequences that could arise from
                            input monomers within the constraints set]
    """

    # copy monomers list for combinatorial addition of monomers to elongating
    # sequences
    sequences = copy.deepcopy(monomers)

    # start building sequences of length i+1 for whole range of i
    for i in range(0,max_length-1):

        # build new sequences of length i+1 by adding monomers on to pre-made
        # sequences of length i
        for monomer in monomers:
            sequences.extend(
                [seq+monomer for seq in sequences
                if len(seq) == i +1])

    # check whether sequencing or just compositional dict is required; if
    # not sequencing, return compositions
    if not sequencing:
        sequences = ["".join(sorted(seq)) for seq in sequences]

    # remove duplicate sequences and / or compositions
    sequences = list(set(sequences))
    
    # remove sequences and / or compositions that do not exceed minimum
    # sequence length threshold specified
    sequences = [seq for seq in sequences if len(seq) >= min_length]
    
    # check for chain terminating monomers; if present, remove sequences that
    # can only be produced by elongating past position of terminating monomers
    if chain_terminators:
        sequences = [
            seq for seq in sequences
            if (seq[0] in chain_terminators
            or seq[-1] in chain_terminators)
        ]

    # check for 0 terminal tags to add to sequences; add if specified
    if start_tags:
        tagged_seqs = []
        for start_tag in start_tags:
            if start_tag in MONOMERS:
                tagged_seqs.extend(
                    [start_tag + sequence for sequence in sequences]
                )
        sequences = tagged_seqs

    # check for -1 terminal tags to add to sequences; add if specified
    if end_tags:
        tagged_seqs = []
        for end_tag in end_tags:
            if end_tag in MONOMERS:
                tagged_seqs.extend(
                    [sequence + end_tag for sequence in sequences]
                )

        sequences = tagged_seqs

    # again, check whether to return sequences or just compositions
    if not sequencing:
        sequences = list(set(["".join(sorted(seq)) for seq in sequences]))

    # check for isobaric target sequences; if specified, return only sequences
    # that are isobaric to one or more of those targets 
    if isobaric_targets:
        isobaric_targets = [
            sorted(target) for target in isobaric_targets
        ]

        sequences = [
            seq for seq in sequences if sorted(seq) in isobaric_targets
        ]
    # return list of sequences and / or compositions
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
                        mass for mass in final_masses]

                    # return sorted list of loss products, removing duplicates
                    return sorted(
                        list(
                            set(
                                [float(f'{mass:.4f}') for mass in final_masses]
                                )
                            )
                        )

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

    # get list of unique monomers within sequence; if only one kind of monomer
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
    sorted_sequences = list(
        set(
            ["".join(sorted(sequence)) for sequence in sequences]
        )
    )
    print(f'organising {len(sequences)} sequences into {len(sorted_sequences)} isobaric groups')

    # generate dictionary of sequence compositions (key = sorted sequences)
    # and lists of sequences matching those compositions (value = sequence list)
    isobaric_dict = {
        sorted_seq: [
            seq for seq in sequences
            if "".join(sorted(seq)) == sorted_seq]
            for sorted_seq in sorted_sequences
        }

    return isobaric_dict

def add_peak_lists_massdict(massdict):
    """
    This function takes an MS1, MS2 or full MSn mass dictionary and generates
    a list of all MS1 and MS2 ions for each sequence in the massdict, adding
    ion masses to a list ("peak_list") which is then used for screening
    spectra.

    Args:
        massdict (dict): dictionary of sequences and corresponding
                subdictionaries of MS1 and /or MS2 ions

    Returns:
        output_dict: same as input massdict, but with extra key-value pair
                added to sequence subdictionaries ({"peak_list": [masses]}
                where masses = list of all MS1 and MS2 ions - including
                monomer-specific signature ions - associated with sequence)
    """
    # initiate output dict to store dictionary to be returned
    output_dict = {}

    # iterate through sequences  and subdicts in input massdict
    for sequence, subdict in massdict.items():

        # initiate list to store all MS1 and MS2 ions associated with sequence
        peak_list = []

        # check sequence subdict for MS1 ions; if present, add to peak_list
        if "MS1" in subdict:
            peak_list.extend(subdict["MS1"])

        # check sequence subdict for MS2 ions; if present, add to peak_list
        if "MS2" in subdict:
            frag_dict = {
                frag: masses
                for frag, masses in subdict["MS2"].items()
                if frag != "signatures" and frag != "unique_fragments"
            }
            for masses in frag_dict.values():
                peak_list.extend(masses)

            # check MS2 subdict for monomer-specific signature ions; if present,
            # add to peak_list
            if "signatures" in subdict["MS2"]:
                signature_fragdict = {
                    frag : masses
                    for frag, masses in subdict["MS2"]["signatures"].items()
                    if frag != "unique_fragments"
                }
                for masses in signature_fragdict.values():
                    peak_list.extend(masses)

        # remove any duplicate m/z values from peak_list
        peak_list = sorted(list(set(peak_list)))

        # add sequence to output dict, with same info as input massdict
        output_dict[sequence] = massdict[sequence]

        # add peak_list subdict to output sequence subdict
        output_dict[sequence].update({"peak_list": peak_list})

    # return massdict with peak lists added for each sequence
    return output_dict

def add_modification_sequence_mass_list(
    mass_list,
    modification_mass,
    mod_mass_diff,
    universal_shift
):
    """
    Takes a mass (or list of masses) corresponding to a sequence and / or
    sequence fragment and adds the mass of a modification on to the sequence
    mass(es), returning a list of modified masses

    Args:
        mass_list (float or list of floats): mass(es) of unmodified sequences
                    and / or sequence fragments
        modification_mass (float): mass of modification to add to unmodified
                    mass(es)
        mod_mass_diff (float): mass lost upon addition of modification to
                    sequence and / or fragments (e.g. water is lost when fatty
                    acids acylate peptides)
        universal_shift (bool): specifies whether modification shifts ALL
                    unmodified masses; if True, only list of modified masses
                    is returned; if False, returned list includes both
                    modified and unmodified masses

    Returns:
        modified_masses (list): list of masses and / or unmodified masses
    """

    # check if mass_list is a list or single mass; if single mass, make list
    if type(mass_list) != list:
        mass_list = [mass_list]

    # make list of modified masses
    modified_masses = [
        round(mass + modification_mass - mod_mass_diff, 4)
        for mass in mass_list
    ]

    # if every mass is to be shifted by the modification, return list of only
    # modified masses, removing duplicates
    if universal_shift:
        return list(set(modified_masses))

    # if universal_shift=False, return list of unmodified masses + modified
    # masses
    modified_masses.extend(mass_list)

    return list(
        set(
            [round(mass, 4) for mass in modified_masses]
        )
    )

def add_terminal_modification_sequence_string(
    sequence,
    modification,
    terminus
):
    """
    Adds three letter code for a terminal modification to a standard sequence
    string

    Args:
        sequence (str): sequence string comprised of monomer one letter codes
        modification (str): three letter code for terminal modification
        terminus (int): specifies terminus to add the sequence; either 0 or -1
        for start terminus and end terminus, respectively

    Returns:
        str : modified sequence string
    """
    if terminus == 0:
        return f"{modification}={sequence}"
    elif terminus == -1:
        return f"{sequence}={modification}"
