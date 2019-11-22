"""
This file contains functions which are essential to all
in silico operations
"""
from ..Constants.GlobalChemicalConstants import *
from ...parameter_handlers import *
import copy
import time
import itertools

# redefine config file variables as their corresponding floats
MASS_DIFF = FG[MASS_DIFF]

def find_sequence_mass(
    sequence,
    mod_markers={
        '(' : ')'
        ,
        'terminal': ['[', ']']
    },
    ms_level=1,
    n_rounded=4
):
    """
    Takes a sequence string and returns neutral monoisotopic mass for sequence

    Arguments:
        sequence {str} -- sequence string made up of constituent monomer one
            letter codes
        mod_markers (dict, optional): dict specifying characters that initiate
            sidechain and terminal modifications. 'sidechain' subdict keys
            mark start of sidechain modification substrings, 'sidechain' subdict
            values mark end of sidechain modifications. 'terminal' subdict
            value = two element list, with element 0 marking start of terminal
            modification substring and element 1 marking end of terminal
            modification. Defaults to {
                'sidechain': {
                    '(': ')
                },
                'terminal': ['[', ']']
            }
        ms_level (int, optional): specifies ms level for calculating. NOTE:
            this will be irrelevant for the vast majority of cases, but has 
            been included in case of future edge cases. Defaults to 1. 
        n_rounded (int, optional): specified number of decimal places to round
            output sequence mass to. Defaults to 4. 

    Returns:
        float: neutral monoisotopic mass of input sequence
    """
    # retrieve list of monomers
    monomers = return_monomers_sequence(
        sequence=sequence,
        mod_markers=mod_markers,
        return_modified=True,
        return_set=False
    )

    # retrieve any terminal modifications from sequence string
    terminal_modifications = remove_substrings_sequence(
        sequence=sequence,
        substring_markers={
            mod_markers['terminal'][0]: mod_markers['terminal'][1]
        },
        return_substrings=True
    )

    # set mass of any terminal modifications to 0; this will be updated 
    # if terminal modifications are found
    terminal_mod_mass = 0

    # check whether terminal modifications have been found; if so, work out
    # neutral mass of all terminal modifications in input sequence 
    if terminal_modifications: 
        
        # sum the neutral masses of all terminal modifications in sequence
        terminal_mod_mass = sum(
            [
                MODIFICATIONS[mod]['mass']
                for mod in terminal_modifications
            ]
        )

        # substract the appropriate ms_level mass diffs for each modification
        # to get the final mass of terminal modifications in the sequence
        terminal_mod_mass -= sum(
            [
                MODIFICATIONS[mod]['mass_diff'][f'ms{ms_level}']
                for mod in terminal_modifications
            ]
        )
    
    # init sequence mass as 0 
    sequence_mass = 0

    for i in range(0, len(monomers)):

        if i == 0: 
            subtract_massdiff = False 
        else:
            subtract_massdiff = True 
        
        sequence_mass += find_monomer_mass(
            monomer=monomers[i],
            subtract_massdiff=subtract_massdiff
        )
    
    # return sequence mass, rounded to specified number of decimal places
    return round(sequence_mass, int(n_rounded))

def find_monomer_mass(
    monomer,
    mod_markers={
        '(': ')',
        'trim_chars': ['~[', ']']
    },
    subtract_massdiff=True,
    ms_level=1
):

    if monomer[0] in ['(', ')']:
        raise Exception(f'{monomer}')
    # if monomer string is only 1 character, no mods are present - return
    # unmodified monomer mass 
    if len(monomer) == 1:

        # retrieve full neutral monoisotopic mass for monomer from polymer-
        # specific config file
        monomer_mass = MONOMERS[monomer][0]

        # check whether to subtract massdiff for monomers
        if subtract_massdiff:
            monomer_mass -= MASS_DIFF
        
        return monomer_mass

    # remove extra characters from monomer string that do not directly 
    # affect monomer mass (e.g. substrings to denote crosslinks with other
    # monomers in full sequence string)
    monomer = "".join([
        c for i, c in enumerate(monomer)
        if c not in mod_markers['trim_chars']
    ])

    # get list of monomer modifications from input monomer string
    monomer_mods = [
        monomer[monomer.find(key)+1:monomer.find(value)]
        for key, value in mod_markers.items()
        if key != 'trim_chars'
    ]
    
    # work out neutral mass of monomer sidechain modifications by summing their
    # neutral masses and subtracting mass_diffs. Info on modifications is 
    # defined in MODIFICATIONS dict found in polymer-specific config file
    monomer_mod_mass = 0
    for mod in monomer_mods: 

        mod_mass = float(MODIFICATIONS[mod]['mass'])
        mod_mass_diff =  MODIFICATIONS[mod]['mass_diff'][f'ms{ms_level}']
        
        try: 
            mod_mass_diff = float(mod_mass_diff)
        except ValueError:
            if mod_mass_diff[0] == '-':
                mod_mass_diff = -FG[mod_mass_diff[1::]]
            else:
                mod_mass_diff = FG[mod_mass_diff]
        
        mod_mass -= mod_mass_diff
        
        monomer_mod_mass += mod_mass
        

    # work out final monomer mass by adding mass of unmodified monomer (found
    # by retrieving monomer one-letter code from first char in monomer string)
    # and mass of modifications
    monomer_mass = float(MONOMERS[monomer[0]][0]) + monomer_mod_mass
    
    
    return monomer_mass

def return_sequence_terminal_modifications(
    sequence
):

    terminal_mods = {}

    if sequence[0] == '[':

        terminal_mods["0"] = sequence[1:sequence.find(']')]
    
    if sequence[-1] == ']':

        mod_start = max([i for i,c in enumerate(sequence) if c == '['])
        terminal_mods["-1"] = sequence[mod_start+1:len(sequence)-1]
    
    return terminal_mods 
    
def remove_substrings_sequence(
    sequence,
    substring_markers={
        '(': ')',
        '~[': ']',
        '[': ']'
    },
    return_substrings=False
):
    """
    Takes a sequence and removes substrings corresponding to modifications. 
    
    Args:
        sequence (str): sequence string comprised of monomers one letter codes
            and / or modification substrings 
        substring_markers (dict, optional): dictionary of chars marking the start of 
            substrings and chars marking their end. Defaults to {
                '(': ')', 
                '~[': ']',
                '[' : ']'
            }
        return_substrings (bool, optional): specify whether to return substrings
            or sequence with substrings removed. if False, sequence is returned;
            if True, list of substrings are returned
    
    Returns:
        str: sequence string with substrings removed
    """
    # init list to store indices in string that are within modificiation
    # substrings
    mod_indices = []

    # init list to store substrings
    substrings = []

    # iterate through sequence, identifying modification substrings and adding
    # their indices to mod_indices
    for i, c in enumerate(sequence):

        if c in substring_markers: 
            
            mod_end = min([
                j for j, d in enumerate(sequence)
                if d == substring_markers[c] and j > i
            ])

            mod_indices.extend([x for x in range(i, mod_end+1)])
            substring = sequence[i:mod_indices[-1]+1]
            substrings.append(substring)

    # check whether substrings are to be returned
    if return_substrings:
        return substrings

    # create new, trimmed sequence from indices that are NOT within modification
    # substrings
    sequence = "".join([
        c for i, c in enumerate(sequence)
        if i not in mod_indices
    ])

    return sequence 

def return_monomers_sequence(
    sequence, 
    mod_markers={
        '(': ')',
        '~[': ']',
        'terminal': ['[', ']']
    },
    return_modified=True,
    return_set=True
):
    """
    Takes a sequence comprised of monomer one-letter codes + / - side chain
    and terminal modification strings, and returns list of monomers in the
    sequence. NOTE: monomers returned can either be exclusively one letter codes
    or include any sidechain modifications (see Args: return_modified).
    
    Args:
        sequence (str): sequence string comprised of monomer one-letter codes
            and / or sidechain and terminal modifications
        mod_markers (dict, optional): dict of characters marking the start 
            of modifications and their corresponding termination characters. 
            All substrings inbetween these characters are assumed to correspond
            to modifications. NOTE: terminal modifications are stored in format:
            {'terminal': [x, y]} where x, y = characters opening and closing
            terminal modification string, respectively. 
             Defaults to {
                 '(': ')',
                 '~[': ']',
                 'terminal': ['[', ']']
            }.
        return_modified (bool, optional): specifies whether to return monomers
            with associated sidechain modification strings. Defaults to True.
        return_set (bool, optional): specifies whether to return list of 
            unique monomers (i.e. list of a set) or not. Defaults to True.
    
    Returns:
        list: list of monomers within input sequence string. 
    """
    # check whether to include sidechain modifications in returned monomer 
    # list. if not, return only 'core' unmodified monomers 
    if not return_modified:
        core_sequence = return_core_sequence(
            sequence=sequence,
            mod_markers=mod_markers
        )
        if return_set:
            return list(set([c for i,c in enumerate(core_sequence)]))
        return [c for i,c in enumerate(core_sequence)]

    # init list to store monomers
    monomers = []
    
    # init list to store indices in sequence which are part of modification
    # substrings
    mod_indices = []

    try:
        # check whether sequence has 0 terminal modification; if so, remove 
        if sequence[0] == mod_markers['terminal'][0]: 

            mod_end = min([
                i for i,c in enumerate(sequence)
                if c == mod_markers['terminal'][1]
            ])

            sequence = sequence[mod_end+1::]
        
        # check whether sequence has -1 terminal modification; if so, remove
        if sequence[-1] == mod_markers['terminal'][-1]:

            mod_start = max(
                i for i,c in enumerate(sequence)
                if c == mod_markers['terminal'][0]
            )

            sequence = sequence[0:mod_start]
    except KeyError:
        pass
    # iterate through sequence, retrieving sidechain-modified monomers and 
    # adding them to list of monomers
    for i, c in enumerate(sequence):

        if c in mod_markers: 
            
            mod_end = min([
                j for j, d in enumerate(sequence)
                if d == mod_markers[c] and j > i
            ])

            # append monomer and its index to monomers list in format [m, i]
            # where m = monomer string and i = index of its first character in
            # sequence
            monomers.append([sequence[i-1:mod_end+1], i-1])
            
            # get indices of all characters containing modified monomer 
            # substrings so that unmodified monomers can be identified
            mod_indices.extend([
                x for x in range(i-1, mod_end+1)
            ])

    # extend monomers list with non-sidechain modified monomers and the index
    # positions in sequence    
    monomers.extend([
        [c, i] for i, c in enumerate(sequence)
        if i not in mod_indices])

    # sort monomers by order they appear in sequence string and remove index
    # markers
    monomers = sorted(
        monomers,
        key=lambda x: x[1]
    )
    monomers = [x[0] for x in monomers]
    
    # return list of unique monomers if return_set is set to True
    if return_set:
        return list(set(monomers))
    
    return monomers

def generate_monomers_list_with_sidechain_mods(
    standard_monomers,
    sidechain_mods,
    universal_modification
):
    """
    This function takes a list of standard, unmodified monomers (i.e. one-letter
    codes), dict of covalent modifications to be added to one or more of the
    monomers, and returns a list of modified and / or unmodified monomers. 
    
    Args:
        standard_monomers (list): list of monomer one letter codes
        sidechain_mods (dict): dict of monomers targeted for modification and
            the three letter codes of modifying groups. Example: 
            {
                'D' : ['Bnz']
            } 
            targets monomer 'D' for modification with 'Bnz'. NOTE: 'Bnz' must 
            be in MODIFICATIONS in polymer-specific config file. 
        universal_modification (bool): specifies whether every monomer targeted
            for modification at sidechain MUST be modified. True if this is
            the case, otherwise False.
    
    Returns:
        list: list of modified and unmodified monomers 
    """

    # return standard input monomers if no sidechain modifications are 
    # specified
    if not sidechain_mods:
        return list(set(standard_monomers))
    
    # init list to store modified monomers
    mod_monomers = []

    # iterate through sidechain modifications dict, adding modification(s)
    # to each target monomer
    for target_monomer, mods in sidechain_mods.items():

        # iterate through each modification, work out whether it's compatible
        # with target monomer; if so, add it to mod_monomers
        for mod in mods:
            compatible_monomers = MODIFICATIONS[mod]["side chain attachments"] 
            if target_monomer not in compatible_monomers:
                print(f'{target_monomer} is not compatible with {mod}')
                print(f'please either update MODIFICATIONS in config file or')
                print(f'choose another monomer-side chain modification')

            else:
                mod_monomers.append(f'{target_monomer}({mod})')
    
    # check whether target monomers are universally modified; if so, remove
    # unmodified monomers that have been targeted 
    if universal_modification: 

        standard_monomers = [
            monomer for monomer in standard_monomers
            if monomer not in [mod_mon[0] for mod_mon in mod_monomers]
        ]

    # add remaining unmodified monomers to mod_monomers, combining all 
    # monomers in one list
    mod_monomers.extend(standard_monomers)
    
    # return list of monomers, both modified and unmodified 
    return (list(set(mod_monomers)))

def return_core_sequence(
    sequence,
    mod_markers={
        '(': ')',
        '~[': ']',
        'terminal': ['[', ']']
    }
):
    """
    Takes a sequence modified with sidechain and/or terminal modifications
    and returns core sequence, with ALL extra characters corresponding to 
    modifications removed. Example: '[Ole]S(Trt)AA' -> 'SAA'
    
    Args:
        sequence (str): sequence string comprised of monomer one letter codes 
            and / or modification strings 
        mod_markers (dict, optional): dictionary of characters marking start
            and end of modification substrings, including terminal modifications. #
            Defaults to {'(': ')','~[': ']','terminal': ['[', ']']}.
    
    Returns:
        str: sequence string comprised solely of monomer one-letter codes, with
            modification substrings removed. 
    """

    # init core sequence string
    core_sequence = ''

    # check whether sequence opens with a 0 terminal modification, remove
    if sequence[0] == mod_markers['terminal'][0]: 

        mod_end = min(
            [
                j for j, d in enumerate(sequence)
                if d == mod_markers['terminal'][1]
            ]
        )
        sequence = sequence[mod_end+1::]
    
    # check whether sequence ends with a -1 terminal modification, remove 
    if sequence[-1] == mod_markers['terminal'][1]:

        mod_start = max(
            [
                j for j,d in enumerate(sequence)
                if d == mod_markers['terminal'][1]
            ]
        )
        sequence = sequence[0:mod_start]

    mod_indices = []
    
    for i, c in enumerate(sequence):

        if c in mod_markers: 

            mod_end = min(
                [
                    j for j,d in enumerate(sequence)
                    if d == mod_markers[c] and j > i
                ]
            )
            mod_indices.extend([j for j in range(i, mod_end+1)])
    
    core_sequence = "".join(
        [
            c for i,c in enumerate(sequence)
            if i not in mod_indices
        ]
    )

    return core_sequence

def add_adducts_sequence_mass(
    mass,
    adducts,
    min_z=1,
    max_z=None,
    mode='pos',
    incoming_charge=0
):
    """
    This functions adds charged adducts to a neutral sequence mass.

    Arguments:
        mass (float) -- neutral monoisotopic mass of sequence OR m/z if 
            incoming_charge > 0.
        adducts (list) -- list of adduct strings. All adducts must be found
                        in either ANIONS or CATIONS dicts in
                        GlobalChemicalConstants.py.
    Keyword Arguments:
        min_z (int) -- minimum ABSOLUTE charge of adduct (default: {1})
        max_z (int) -- maximum ABSOLUTE charge of adduct. If this is set to
                        None, maximum charge is assumed to be maximum oxidation
                        state of adduct ions (default: {None}).
        mode (str) -- either 'pos' or 'neg' for positive and negative mode
                        mass spectrometry, respectively (default: {'pos'}).
        incoming_charge (int, optional): specifies charge associated with 
            input mass. This will only be greater than 0 for MS2 fragments 
            with intrinsic charges (e.g. peptide b fragments). Defaults to 0.
    Returns:
        charged_sequence_masses (list) -- list of m/z values for charged sequences
                        with adducts.
    """
    # list anions and cations within adducts
    anions = [adduct for adduct in adducts if adduct in ANIONS]
    cations = [adduct for adduct in adducts if adduct in CATIONS]

    if type(mass) != list:
        masses = [mass]
    else:
        masses = mass

    # check whether counterions need to be considered for adducts which have the
    # opposite charge from mode. If so, return charged adducts from adduct_complex
    # function
    if (mode == 'pos' and len(anions) > 0) or (mode == 'neg' and len(cations) > 0):
        charged_sequence_masses = add_adduct_complexes_sequence_mass(
            sequence='',
            neutral_mass=mass,
            adducts=adducts,
            min_z=min_z,
            max_z=max_z)
        return charged_sequence_masses

    # retrieve adduct masses and charges from GlobalChemicalConstants
    if mode == 'pos':
        ions = CATIONS
    elif mode == 'neg':
        ions = ANIONS

    # initiate list to add charged adduct m/z values
    charged_sequence_masses = []

    # iterate through adducts and return m/z values of charged species
    # include multiply charged species that fit within constrains of min_z, max_z
    # and ion oxidation states given in GlobalChemicalConstants
    
    for adduct in adducts:
        
        adduct_mass = ions[adduct][0]
        min_charge = int(ions[adduct][1])
        max_charge = int(ions[adduct][2])
    
            
        # get maximum and minimum charge - either specified or use default
        # minimum and maximum charge states of adducts from
        # GlobalChemicalConstants
        
        if not max_z:
            maximum_charge = max_charge
        else:
            maximum_charge = min([int(max_z), max_charge])

        # add more comments
        min_z = max([min_z, min_charge])

        for i in range(min_z, maximum_charge+1):
            
            charged_sequence_masses.extend(
                [(mass+adduct_mass)/i for mass in masses])

    
    #finally, return list of charged sequence m/z values
    return sorted(list(set(
                [float(f"{mass:0.4f}")
                for mass in charged_sequence_masses])))

def find_functional_groups_monomer(
    monomer
):
    """ This function returns the functional groups present in a monomer.

    Arguments:
        monomer (string) -- One letter monomer code that has associated neutral
        monoisotopic mass, list of functional groups with their associated
        frequency per monomer.

    Returns:
        func_groups -- List of functional groups present in the monomer.
    """

    monomer_info = MONOMERS[monomer]

    func_groups = [x[0] for x in monomer_info[1]]

    return func_groups

def reverse_sequence(
    sequence,
    mod_markers={
        '(': ')',
        'terminal': ['[', ']']
    }
):
    """
    Takes a sequence comrpised of monomer one-letter codes and / or any 
    terminal and sidechain modification substrings and outputs the reverse 
    of that sequence. 
    Example: 
        '[Ole]S(Trt)GAVS(Trt)' -> 'S(Trt)VAGS(Trt)[Ole]'
    
    Args:
        sequence (str): sequence string comprised of monomer one letter codes 
            and / or modification substrings 
        mod_markers (dict, optional): dict to specify which characters mark
            beginning and end of sidechain and terminal modification substrings. 
            Defaults to {'(': ')','terminal': ['[', ']']}.
    
    Returns:
        str: reverse of input sequence
    """
    reversed_sequence = sequence[::-1]

    # reverse NON-TERMINAL mod markers to deal with reversed sequence
    r_mod_markers = {
        v : k
        for k, v in mod_markers.items()
        if k != 'terminal'
    }

    # if terminal mod markers have been specified in mod_markers, reverse these
    # also
    try:
        terminal_mod_markers = mod_markers['terminal']
        r_terminal_markers = {
            terminal_mod_markers[1]: terminal_mod_markers[0]
        }
    except KeyError:
        r_terminal_markers = {}

    # generate list of reversed modification strings from reversed sequence
    # this list will only contain reversed NON-TERMINAL mods
    reversed_mods = remove_substrings_sequence(
        sequence=reversed_sequence,
        substring_markers=r_mod_markers,
        return_substrings=True
    )

    # generate list of reversed terminal mods from reversed sequence
    reversed_terminal_mods = remove_substrings_sequence(
        sequence=reversed_sequence,
        substring_markers=r_terminal_markers,
        return_substrings=True
    )
    
    mod_monomers = {}
    
    # check whether any reversed mods have been found in reversed sequence
    # if so, switch them back to regular non-reversed substrings
    if reversed_mods: 
        
        # iterate through each reversed side chain modification and ensure 
        # the now-reversed modification substring is replaced with the original
        # non-reversed string
        for r_mod in reversed_mods:
            reversed_sequence = reversed_sequence.replace(r_mod, r_mod[::-1])
            
            # find positons of side chain modification substrings in reversed
            # sequence along with their associated monomers
            for i, c in enumerate(sequence):
                if c == r_mod[::-1][0]:
                    
                    if sequence[i:i+len(r_mod)] == r_mod[::-1]:
                        r_mod_monomer = sequence[i-1:i+len(r_mod)]
                        non_r_mod_monomer = f'{r_mod_monomer[1::]}{r_mod_monomer[0]}'
                        mod_monomers[r_mod_monomer] = non_r_mod_monomer

    # move position of sidechain modification substrings in reversed sequence
    # to ensure they remain associated with correct monomers
    if mod_monomers: 
        for monomer, non_r_monomer in mod_monomers.items(): 
            reversed_sequence = reversed_sequence.replace(non_r_monomer, monomer)
   
    # check whether any terminal mods have been found in reversed sequence
    # if so, switch them back to regular non-reversed substrings
    if reversed_terminal_mods: 

        for r_t_mod in reversed_terminal_mods:
            reversed_sequence = reversed_sequence.replace(r_t_mod, r_t_mod[::-1])
    
    return reversed_sequence

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
    metal centres, ions with counterions and associated solvent species etc.

    Args:
        sequence (str): sequence string comprised of monomer one letter codes.
        neutral_mass (float): neutral monoisotopic mass of sequence.
        adducts (list): list of adduct strings.
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

    Returns:
        masses (list): list of masses of sequences with adducts.
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
    start_tags=None,
    end_tags=None,
    isobaric_targets=None):
    """
    This function takes a list of input monomers and outputs all possible
    sequences or compositions that could arise within the constraints set,
    which are described below.

    Arguments:
        monomers (list) -- list of monomer one letter codes.
        max_length (int) -- maximum sequence length (in monomer units).
    Keyword Arguments:
        min_length (int) -- minimum sequence length (in monomer units).
                            (default: {1}).
        sequencing (bool) -- specifies whether all possible sequences are to
                            be enumerated or just compositions. Set to False if
                            you just need to screen for compositions]
                            (default: {True})
        start_tags {list} -- list of monomers. If this list is input, sequences
                            are tagged at terminus 0 by each of the monomers
                            in start_tags - one tag per sequence, and only
                            tagged sequences are returned. (default: {None}).
        end_tags {list} -- list of monomers. If this list is input, sequences
                            are tagged at terminus -1 by each of the monomers
                            in end_tags - one tag per sequence, and only
                            tagged sequences are returned. (default: {None}).
        isobaric_targets {list} -- list of sequences and / or compositions that
                            output sequences must be isobaric to. (default: 
                            {None}).
    Returns:
        sequences {list} -- list of possible sequences that could arise from
                            input monomers within the constraints set.
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
                [
                    seq+monomer for seq in sequences
                    if (len(return_monomers_sequence(
                        sequence=seq,
                        return_modified=True,
                        return_set=False)) == i +1)
                ]
        )
    
    # check whether unique sequences or just compositions are to be returned 
    if not sequencing: 
        sequences = [
            "".join(sorted(return_monomers_sequence(
                sequence=sequence,
                return_modified=True,
                return_set=False)))
                for sequence in sequences
        ]
    
    # remove duplicate sequences and / or compositions
    sequences = list(set(sequences))
    
    # remove sequences and / or compositions that do not exceed minimum
    # sequence length threshold specified
    sequences = [seq for seq in sequences if len(seq) >= min_length]
    
    # check for 0 terminal tags to add to sequences; add if specified
    if start_tags:
        tagged_seqs = []
        for start_tag in start_tags:
            if start_tag in MONOMERS:
                tagged_seqs.extend([
                    start_tag + sequence for sequence in sequences])
                
        sequences = tagged_seqs
    
    # check for -1 terminal tags to add to sequences; add if specified
    if end_tags:
        tagged_seqs = []
        for end_tag in end_tags:
            if end_tag in MONOMERS:
                tagged_seqs.extend(
                    [sequence + end_tag for sequence in sequences])

        sequences = tagged_seqs
    
    # again, check whether to return sequences or just compositions
    if not sequencing:
        sequences = [
            "".join(sorted(return_monomers_sequence(
                sequence=sequence,
                return_modified=True,
                return_set=False)))
                for sequence in sequences
        ]

    
    # check for isobaric target sequences; if specified, return only sequences
    # that are isobaric to one or more of those targets 
    if isobaric_targets:
        isobaric_targets = [
            sorted(target) for target in isobaric_targets]

        sequences = [
            seq for seq in sequences if sorted(seq) in isobaric_targets]
    
    # return list of sequences and / or compositions
    return sequences


def generate_all_sequences_rxn_classes(monomers):
    """
    This function will be written to build sequence lists in cases where
    monomers are not universally cross-reactive - i.e. more complicated
    operations are required.

    THIS FUNCTION IS NOT FINISHED

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
    Takes a sequence, mass and subtracts monomer-specific side chain neutral
    loss products as specified in LOSS_PRODUCTS dictionary found in polymer
    config file.

    Arguments:
        sequence (str) -- sequence string comprised of monomer one letter
                            codes.
        sequence_masses (float or list) -- float for a single mass, or list
                                    of floats for multiple associated masses,
                                    e.g. for various adducts.

    Keyword Arguments:
        max_total_losses (int) -- specifies the maximum number of loss
                                products to be subtracted from a given sequence
                                if None, all possible loss products will be
                                incorporated into final mass list.
                                (default: {None}).

    Returns:
        final_masses (list) -- list of m/z values for sequence +/- side
                                chain neutral losses.
    """
    # make sequence_mass a list if not already
    if type(sequence_masses) != list:
        sequence_masses = [sequence_masses]
    
    if max_total_losses == 0:
        return None

    # get dictionary of loss-product prone monomers and associated loss products
    monomers = return_monomers_sequence(
        sequence=sequence,
        return_modified=True,
        return_set=True
    )

    # get dict of loss product-prone monomers and their losses 
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

            # check if loss is a string, if so redefine as the corresponding float
            if type(loss)==str:
                loss = FG[loss]
            for i in range(1, occurence+1):
                neutral_loss = loss*i
                monomer_losses.extend([
                    mass-neutral_loss for mass in final_masses])

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
                    return sorted(list(set(
                        [round(float(mass), 4) for mass in final_masses])))
                               
        # add monomer-specific loss products to final sequence masses,
        # remove duplicate final masses
        # make masses four point floats
        final_masses.extend(monomer_losses)
        final_masses = sorted(list(set(
            [round(float(mass), 4) for mass in final_masses])))

    return final_masses

def generate_reading_frames_sequence(
    sequence,
    mod_markers={
        '(':')',
        '~[': ']'}
    ):
    """
    Takes a linear sequence and outputs a list of reading frame shifts for
    later use in generating proposed fragments for cyclic sequences. A reading
    frame shift is carried out by taking the final monomer in the sequence
    and making it the first monomer, shifting the remaining monomers up 1 index
    in the polymer chain. Example: one reading frame shift of sequence 'AAV'
    produces shifted sequence 'VAA'.

    Args:
        sequence (str): sequence string consisting of monomer one letter codes.

    Returns:
        reading_frames: list of unique reading frames.
    """
    
    # get list of monomers, sorted by order they appear in sequence
    monomers = return_monomers_sequence(
        sequence=sequence,
        mod_markers=mod_markers,
        return_modified=True,
        return_set=False
    )

    # return no reading frames if there is only one type of monomer
    if sorted(list(set(monomers))) == sorted(monomers):
        return []
    
    # init list to store unique reading frames 
    reading_frames = [sequence]

    for i in range(0, len(monomers)):
        
        # generate new reading frame by making last monomer first monomer
        # in sequence
        r_frame = f"{''.join(monomers[0:len(monomers)-1])}{monomers[-1]}"

        monomers = return_monomers_sequence(
            sequence=r_frame,
            mod_markers=mod_markers,
            return_modified=True,
            return_set=False
        )

        # add reading frame to list of reading frames if it is unique
        if r_frame not in reading_frames:
            reading_frames.append(r_frame)

    return reading_frames

def generate_dict_isobaric_sequences(sequences):
    """
    Takes a list of sequences and groups isobaric sequences in a dictionary.

    Args:
        sequences (list): list of sequence strings

    Returns:
        dict -- dictionary of isobaric groups (key = sorted sequence, value = list
            of sequences isobaric to sorted sequence).
    """
    # get list of sequence compositions
    sorted_sequences = list(set(
        ["".join(sorted(sequence)) for sequence in sequences]))

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
                subdictionaries of MS1 and /or MS2 ions.

    Returns:
        output_dict: same as input massdict, but with extra key-value pair
                added to sequence subdictionaries ({"peak_list": [masses]}
                where masses = list of all MS1 and MS2 ions - including
                monomer-specific signature ions - associated with sequence).
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
                if (frag != "signatures" 
                and frag.find("mod") == -1 
                and frag != "unique_fragments")
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
                for sig, masses in signature_fragdict.items():
                    if sig != 'terminal_modifications':
                        peak_list.extend(masses)
                    else:
                        for t_masses in masses.values():
                            peak_list.extend(t_masses)
                

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
    mass(es), returning a list of modified masses.

    Args:
        mass_list (float or list of floats): mass(es) of unmodified sequences
                    and / or sequence fragments.
        modification_mass (float): mass of modification to add to unmodified
                    mass(es).
        mod_mass_diff (float): mass lost upon addition of modification to
                    sequence and / or fragments (e.g. water is lost when fatty
                    acids acylate peptides).
        universal_shift (bool): specifies whether modification shifts ALL
                    unmodified masses; if True, only list of modified masses
                    is returned; if False, returned list includes both
                    modified and unmodified masses.

    Returns:
        modified_masses (list): list of masses and / or unmodified masses.
    """
    # check if mass_list is a list or single mass; if single mass, make list
    if type(mass_list) != list:
        mass_list = [mass_list]

    try:
        mod_mass_diff = float(mod_mass_diff)
    except ValueError:
        if mod_mass_diff[0] == "-":
            mod_mass_diff = -FG[mod_mass_diff[1::]]
        else:
            mod_mass_diff = FG[mod_mass_diff]
    # make list of modified masses
    modified_masses = [
        round(mass + float(modification_mass) - mod_mass_diff, 4)
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
    terminal_modification,
    terminus
):
    """
    Adds three letter code for a terminal modification to a standard sequence
    string to return a single modified sequence string OR returns a list of
    terminally modified sequences if more than one terminal_modification 
    is given (see Args: terminal_modification)

    Args:
        sequence (str): sequence string comprised of monomer one letter codes
        terminal_modification (str or list): single string corresponding to
            a single terminal modification or a list of terminal modification
            strings
        terminus (int or str): specifies terminus to add the sequence; either 
        0/"0" or -1/"-1" for start terminus and end terminus, respectively

    Returns:
        str or list: modified sequence string or list of modified sequence
            strings 
    """

    if str(terminus) == "0":
        
        modified_sequences = [
            f"[{mod}]{sequence}"
            for mod in terminal_modification
        ]
        
    elif str(terminus) == "-1":        
        modified_sequences = [
            f"{sequence}[{mod}]"
            for mod in terminal_modification
        ]
    
    return list(set(modified_sequences))
