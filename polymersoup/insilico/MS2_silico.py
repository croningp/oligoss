"""
This file contains functions for generating theoretical MS2 data for ALife
Polymer Data
"""
from .Constants.GlobalChemicalConstants import *
from .Config_files.Depsipeptide_config import *
from .silico_helpers.insilico_helpers import *

def build_fragment_series_single_sequence(
    sequence,
    fragment_series,
    adducts,
    mod_markers={
        '(': ')',
        'terminal': ['[', ']']
    },
    universal_mass_shift=False,
    mode='pos',
    losses=True,
    max_total_losses=None,
    loss_product_adducts=False,
    check_exceptions=False
):
    
    # define frag_info (makes nested dictionary navigation more readable)
    frag_info = FRAG_SERIES[fragment_series]

    # get fragment terminus (i.e. end of sequence fragment series starts from)
    terminus = frag_info['terminus']
    
    # get terminal modifications at terminus -1 and terminus 0, if any
    terminal_mods = return_sequence_terminal_modifications(
        sequence=sequence
    )

    # reverse sequence if fragment series is indexed from -1 terminus
    if int(terminus) == -1:
        
        sequence = reverse_sequence(
            sequence=sequence,
            mod_markers=mod_markers
        )
        

        
    # get fragment mass diff (i.e. mass difference between fragments and
    # equivalent sub-sequences)
    mass_diff = frag_info['mass_diff'][mode]
    
    # check whether fragment series is possible in current mode. if not, 
    # no fragments will be returned. Fragment series and their compatible modes
    # are defined in polymer-specific config file
    if not mass_diff:
        print(f'{fragment_series} fragment series does not exist for {mode} mode')
        return {}
    
    # get list of monomers contained in sequence, sorted by order they appear
    # in sequence
    monomers = return_monomers_sequence(
        sequence=sequence,
        mod_markers=mod_markers,
        return_modified=True,
        return_set=False
    )

    # init values to add to start and end terminal fragments
    start_terminal_mass = 0
    end_terminal_mass = 0
    
    # check for terminal modifications, and retrieve their masses 
    if terminal_mods:

        for mod_terminus in terminal_mods: 
            
            print(f'terminal mods, sequence = {terminal_mods}, {sequence}')
            
            # get info for terminal modification
            terminal_mod_info = MODIFICATIONS[terminal_mods[str(mod_terminus)]]

            try:
                mod_mass = float(terminal_mod_info["mass"])
            except ValueError:
                if mod_mass[0] == '-':
                    mod_mass = -FG[mod_mass[1::]]
                else:
                    mod_mass = FG[mod_mass]
            try:
                mod_mass -= float(terminal_mod_info["mass_diff"]["ms2"])
            except ValueError:
                mod_mass_diff = terminal_mod_info["mass_diff"]["ms2"]
                if mod_mass_diff[0] == '-':
                    mod_mass_diff = -FG[mod_mass_diff[1::]]
                else:
                    mod_mass_diff = FG[mod_mass_diff]
                mod_mass -= mod_mass_diff

            if int(mod_terminus) == int(terminus):
                start_terminal_mass += mod_mass 
            else: 
                end_terminal_mass += mod_mass
    
    # get sequence length from number of constituent monomers
    sequence_length = len(monomers)

    # get fragment series starting position along linear polymer (i.e. number
    # of residues from terminus to start building fragment series)
    start = frag_info['start']

    # get fragment series ending position along linear polymer (i.e. number of
    # residues from terminus to stop building fragment series)
    end = frag_info['end']
    end = sequence_length - end

    # init dict to store fragments and their associated masses
    frag_dict = {}
    
    # iterate through fragment series length, creating fragments from each
    # subsequence in series
    for i in range(max([0, start]), end):

        # create subsequence for current fragment
        subsequence = "".join(monomers[0:i+1])

        current_fragment = f"{fragment_series}{i+1}"
        
        # get neutral monoisotopic mass of subsequence 
        neutral_fragment_masses = [find_sequence_mass(
            sequence=subsequence,
            mod_markers=mod_markers,
            ms_level=2)]
        
        frag_mass_exception = None

        if check_exceptions:
            
            # check for exceptions to fragment mass diff
            frag_mass_exception = apply_fragmentation_exceptions(
                fragment_id=current_fragment,
                sub_sequence=subsequence,
                mode=mode
            )
           
        if not frag_mass_exception:
            
            # subtract mass diff for fragment series (e.g. +H for positive peptide
            # y fragments, -OH for positive peptide b fragments)
            frag_massdiff = FRAG_SERIES[fragment_series]["mass_diff"][mode]
        else:
            frag_massdiff = frag_mass_exception
        
        try:
            frag_massdiff = float(frag_massdiff)
        except ValueError:
            if frag_massdiff[0] == '-':
                frag_massdiff = -FG[frag_massdiff[1::]]
            else:
                frag_massdiff = FG[frag_massdiff]
        
        # apply massdiff to all fragments, creating charged fragment species
        neutral_fragment_masses = [
            mass + frag_massdiff
            for mass in neutral_fragment_masses
        ]

        # if no losses are allowed, set maximum total losses to 0 so no 
        # loss products are returned
        if not losses:
            max_total_losses = 0 
        
        # work out sidechain-specific neutral loss products 
        loss_products = add_sidechain_neutral_loss_products_sequence(
            sequence=subsequence,
            sequence_masses=neutral_fragment_masses,
            max_total_losses=max_total_losses)
        
        # check whether to add extra (i.e. non-intrinsic) adducts to loss 
        # products
        if loss_product_adducts:
            neutral_fragment_masses.extend(loss_products)
        
        # add non-intrinsic adducts to fragment masses 
        current_frag_dict = {current_fragment: neutral_fragment_masses}
        current_frag_dict = add_adducts_ms2_fragments(
            fragment_dict=current_frag_dict,
            adducts=adducts,
            mode=mode,
            return_base_fragments=True
        )

        # make sure all fragment masses (loss products + / - adducts)
        # are accounted for
        final_frag_masses = current_frag_dict[current_fragment]
        
        final_frag_masses.extend(loss_products)        
        
        total_terminal_mass = start_terminal_mass 

        # check whether final fragment in series
        if i == end-1:

            if sequence_length - end == 0: 
                total_terminal_mass += end_terminal_mass 
        
        modified_frag_masses = [
            mass + total_terminal_mass
            for mass in final_frag_masses
        ]

        if universal_mass_shift:
            final_frag_masses = modified_frag_masses
        else:
            final_frag_masses.extend(modified_frag_masses)

        # add fragment and neutral masses to frag_dict
        frag_dict[current_fragment] = list(set(
            [round(mass, 4) for mass in final_frag_masses]))

    # finally, return MS2 fragment dictionary
    return frag_dict

def return_modification_signature_ions(
    sequence,
    mode
):

    terminal_modifications = return_sequence_terminal_modifications(
        sequence=sequence
    )

    if not terminal_modifications:
        return {}
    
    mod_signatures = {}

    for mod in terminal_modifications.values():

        mod_signatures.update(
            {
                mod: MODIFICATIONS[mod]["free_mod_fragments"][mode]
            }
        )

    return mod_signatures

def build_multiple_fragment_series_single_sequence(
    sequence,
    fragment_series,
    adducts,
    mod_markers={
        '(': ')',
        'terminal': ['[', ']']
    },
    mode='pos',
    add_signatures=True,
    signatures=None,
    losses=True,
    max_total_losses=None,
    loss_product_adducts=True,
    universal_mass_shift=False
):
    """
    Takes a single sequence, list of fragment series (e.g. ['b', 'y'] for
    standard peptide b- and y- fragments) and outputs a dictionary of MS2
    fragments and corresponding m/z values.

    Args:
        sequence (str): Sequence string comprised of monomer one letter codes.
        fragment_series (list, optional): List of fragment one letter codes
                        strings; if None, all fragments in FRAG_SERIES from
                        polymer config file are included in final fragment
                        dictionary. Defaults to None.
        mode (str, optional): Either 'pos' or 'neg' to specificy positive or
                        negative mode, respectively. Defaults to 'pos'.
        add_signatures (bool, optional): Specifies whether to add
                        monomer-specific MS2 signature ions to fragment_dict.
                        Defaults to True.
        signatures ([type], optional): List of signature ion strings, to be
                        added to fragment_dict if add_signatures=True. If None,
                        all MS2 signature ions in polymer config file
                        MS2_SIGNATURE_IONS dict are added to fragment_dict.
                        Defaults to None.

    Returns:
        fragment_dict: dictionary of MS2 fragments and corresponding m/z values.
    """
    # if no fragment series are specified, all fragment types in FRAG_SERIES
    # dict from polymer config file are to be incorporated into MS2 fragment
    # dict
    if not fragment_series:
        fragment_series = FRAG_SERIES.keys()

    # initiate fragment dictionary to be populated
    fragment_dict = {}

    # set default to not look for exceptions to standard fragmentation rules
    check_exceptions = False

    # iterate through fragment series list and update fragment dictionary for
    # each fragment type
    for series in fragment_series:
        if 'exceptions' in FRAG_SERIES[series]["mass_diff"].keys():
            check_exceptions = True

        fragment_dict.update(
            build_fragment_series_single_sequence(
                sequence=sequence,
                fragment_series=series,
                adducts=adducts,
                mod_markers=mod_markers,
                universal_mass_shift=universal_mass_shift,
                check_exceptions=check_exceptions)
            )

    # checks if monomer-specific MS2 signature ion fragments are to be added
    # if specified, signature ions are then added to fragment dict
    if add_signatures:
        fragment_dict.update({
            'signatures': generate_monomer_ms2_signature_ions_sequence(
                sequence=sequence,
                signatures=signatures)})

    return fragment_dict

def apply_fragmentation_exceptions(
    fragment_id,
    sub_sequence,
    mode
): 
    
    # get monomers, listed in order they appear in subsequence
    monomers = return_monomers_sequence(
        sequence=sub_sequence,
        return_set=False,
        return_modified=True
    )

    # get fragment series and index from fragment id
    fragment_series = fragment_id[0]

    try:
        # get positions of exceptions and their reactivity classes 
        exception_positions = FRAG_SERIES[fragment_series]["mass_diff"]["exceptions"][mode]
    except KeyError:
        return None

    for pos in exception_positions:
        
        monomer = monomers[int(pos)][0]

        rxn_classes = exception_positions[pos].keys()
        
        for rxn_class in rxn_classes:
            if monomer in REACTIVITY_CLASSES[rxn_class][1]:
                try:
                    exception_mass = float(exception_positions[pos][rxn_class])
                except ValueError:
                    exception_key = exception_positions[pos][rxn_class]
                if exception_key[0] == "-":
                    exception_mass = -FG[exception_key[1::]]
                else:
                    exception_mass = FG[exception_key]
                return exception_mass

    return None

def add_adducts_ms2_fragments(
    fragment_dict,
    adducts,
    mode='pos',
    return_base_fragments=True
):
    """
    This function takes a dictionary of MS2 fragments and associated m/z
    values, and adds singly charged adducts.

    Arguments:
        fragment_dict (dict) -- dictionary of fragments and associated m/z
                            values. keys = fragment id strings, values = lists
                            of associated masses. See function:
                            build_fragment_series_single_sequence.
        adducts (list) -- list of adducts to be added to MS2 fragments.

    Keyword Arguments:
        mode {str} -- either 'pos' or 'neg' for positive and negative mode,
                        respectively (default: {'pos'})
        return_base_fragments {bool} -- specific whether to return input
                        fragment masses in adduct list (default: {True})

    Raises:
        Exception: raised if any of the input adducts do not match the
                    detection mode (i.e. anions in positive mode, cations in
                    positive mode). Functions for complex adducts are yet to
                    be written.

    Returns:
        adduct_fragdict -- dictionary of fragment string codes and associated
                        m/z value lists with adducts. ALL IN +1 CHARGE STATE.
    """
    # list fragments, excluding monomer-specific signature ion fragments
    fragments = [key for key in fragment_dict if key != 'signatures']

    # retrieve fragment series (i.e. fragment one letter codes) for all 
    # fragments in fragment_dict
    fragment_series = list(set(
        [fragment[0] for fragment in fragments]))
    
    # initate dictionary to populate with fragments and associated
    # fragment+adduct masses
    fragment_adducts = {}

    # iterate through fragments by type, adding adducts to each fragment subset
    for series in fragment_series:
        current_fragments = {
            fragment: fragment_dict[fragment]
            for fragment in fragment_dict
            if fragment[0] == series
        }

        # retrieve info on fragment series
        frag_info = FRAG_SERIES[series]

        # get 'intrinsic adducts' for fragment series - i.e. adducts that are
        # added to fragments as default by fragment builder and therefore must
        # be removed before addition of non-standard adducts to fragments
        intrinsic_adduct = 0
        if 'intrinsic_adducts' in frag_info:
            intrinsic_adduct = frag_info['intrinsic_adducts'][mode]
            if type(intrinsic_adduct) == str:
                intrinsic_adduct = FG[intrinsic_adduct]

        # some fragment series have intrinsic charges (not adducts) and are
        # not found with extra adducts in a singly charged state - e.g.
        # peptide b fragments
        if 'permissible_adducts' in frag_info:
            adducts = [adduct for adduct in adducts
            if adduct in frag_info['permissible_adducts'][mode]]

        # organise adducts by type - cationic or anionic
        cations = [adduct for adduct in adducts if adduct in CATIONS]
        anions = [adduct for adduct in adducts if adduct in ANIONS]

        # check that all adducts are compatible with overall charge of
        # fragmenting species - i.e. cations for positive mode, anions for
        # negative. If charge state of adducts does not match overall charge
        # state, counterions must be added to adducts in a separate function
        if mode == 'pos' and sorted(cations) != sorted(adducts):
            raise Exception('function for adduct ions and counterions not done')
        elif mode == 'neg' and sorted(anions) != sorted(adducts):
            raise Exception('function for adduct ions and counterions not done')

        # iterate through fragments in series, adding adduct masses
        for fragment, masses in current_fragments.items():
            adduct_masses = []
            for mass in masses:
                mass -= intrinsic_adduct
                adduct_masses.extend(add_adducts_sequence_mass(
                    neutral_mass=mass,
                    adducts=adducts,
                    min_z=1,
                    max_z=1,
                    mode=mode
                )
            )

            # if standard fragments (i.e. with unchanged input masses) are to
            # be returned, include these in output adduct_masses list
            if return_base_fragments:
                adduct_masses.extend(masses)

            adduct_masses = [float(f'{mass:.4f}') for mass in adduct_masses]

            # add fragment and associated adducts to output adduct_fragdict,
            # removing any duplicate masses - which should not be there anyway
            fragment_adducts[fragment] = sorted(list(set(adduct_masses)))

    # finally, return dictionary of fragments and associated adduct mass lists
    return fragment_adducts

def generate_monomer_ms2_signature_ions_sequence(
    sequence,
    signatures=None,
    mode='pos',
    mod_signatures=True
):
    """
    This function generates a dictionary of monomer signature ions for
    a sequence.

    Arguments:
        sequence (str) -- sequence string with constituent monomer one letter
                            codes.
    Keyword Arguements:
        signatures (list) -- list of signature fragment codes - e.g. ['Im']
                            for amino acid immonium fragments. If None, all
                            signature ion types will be returned in signature
                            ion dictionary
                            (default: {None})
        mod_signatures (bool, optional) -- specifies whether to add signature 
            ions for terminal and / or sidechain covalent modifications. Defaults
            to True.
    Returns:
        dict -- dictionary of monomer signature fragments
                            and associated m/z values
    """
    # get list of unique monomers within sequence
    monomers = return_monomers_sequence(
        sequence=sequence,
        return_modified=False
    )

    # get list of signature ion id strings used to access signature ion values
    # from polymer config file
    if not signatures:
        signatures = MS2_SIGNATURE_IONS.keys()

    # initiate output MS2 signature ion dict
    MS2_signature_dict = {}

    # iterate through signature ion types and add to signature ion dictionary
    for signature in signatures:
        if signature != "dominant":
            signature_ions = MS2_SIGNATURE_IONS[signature]

            MS2_signature_dict.update(
                {
                    f'{signature}{monomer}': signature_ions[monomer]
                    for monomer in monomers
                    if monomer in signature_ions
                }
            )

    modifications = return_sequence_terminal_modifications(
        sequence=sequence
    )

    terminal_dict = {}
    if mod_signatures and modifications: 

        try: 
            terminal_dict = {
                mod: MODIFICATIONS[mod]["free_mod_fragments"][mode]
                for mod in modifications.values()
            }

        except KeyError:
            pass
            
    MS2_signature_dict.update(
        {'terminal_modifications': terminal_dict}
    )

    # finally, return dictionary of monomer signature ions and associated
    # lists of m/z values
    return MS2_signature_dict

def check_monomers_for_exceptions(
    monomers,
    fragment_series,
    mode
):
    """
    This function takes a list of monomers, MSn fragment series and works out
    whether any monomers have fragmentation patterns atypical for that 
    fragment series - e.g. hydroxy acids for positive mode peptide y fragments
    
    Args:
        monomers (list): list of monomer strings (can be one letter codes or
            modified monomer strings)
        fragment_series (list): list of fragment series one letter codes (
            e.g. ['b', 'y'] for peptide b and y fragment series)
        mode (str): ither 'pos' or 'neg' for positive and negative mode 
            mass spectrometry, respectively
    
    Returns:
        bool: True if any exceptions to standard fragmentation rules occur, 
            False if not
    """

    # retrieve fragment series that have exceptions to standard rules in 
    # at least some cases
    exception_frags = [
        frag for frag in fragment_series
        if 'exceptions' in FRAG_SERIES[frag]["mass_diff"]
    ]

    # if no exceptions are possible, return False 
    if not exception_frags:
        return False 
    
    # check whether any fragment series have exceptions to standard rules for
    # mass spec mode used
    exception_frags = [
        frag for frag in exception_frags 
        if mode in FRAG_SERIES[frag]["mass_diff"]["exceptions"]
    ]

    # return False if no fragment series have exceptions in current mass 
    # spec mode
    if not exception_frags:
        return False 
    
    # init list to store reactivity classes and exceptions that apply in 
    # current mass spec mode 
    exception_rxn_classes = []

    # update exception rxn classes for every rxn class that has exception 
    # in current mass spec mode
    for frag in exception_frags:
        
        exception_dict = FRAG_SERIES[frag]["mass_diff"]["exceptions"][mode]
        for terminus in exception_dict:
            exception_rxn_classes.extend(exception_dict[terminus].keys())

    # work out if any monomers have functional group(s) that require exceptions
    # and return True if this is the case
    for rxn_class in exception_rxn_classes:
      
        class_monomers = [
            monomer for monomer in monomers 
            if monomer[0] in REACTIVITY_CLASSES[rxn_class][1]
        ]
        if class_monomers:
            return True

    return False 
    
def generate_ms2_mass_dictionary(
    sequences,
    fragment_series,
    adducts,
    mode,
    add_signatures=True,
    signatures=None,
    losses=True,
    max_total_losses=None,
    loss_product_adducts=False,
    uniques=True
):
    """
    Takes a list of sequences and generates full ms2 fragment dictionary.

    Args:
        sequences (str): sequence string comprised of monomer one letter codes.
        fragment_series (list): list of fragment types - must be keys in
                            FRAG_SERIES dict in polymer config file.
        adducts (list): list of adduct strings.
        mode (str): either 'pos' or 'neg' for positive and negative mode,
                    respectively.
        add_signatures (bool, optional): specifies whether to add
                            monomer-specific MS2 signature ions.
                            Defaults to True.
        signatures (list, optional): list of signature ion types to add; if
                            None, all possible signature ions are added.
                            Defaults to None.
        losses (bool, optional): specifies whether to include monomer-specific
                            side chain loss products. Defaults to True.
        max_total_losses (int, optional): specifies whether to cap number of
                            side chain loss products per fragment. If None,
                            all possible loss products are generated.
                             Defaults to None.
        loss_product_adducts (bool, optional): specify whether to add non-
            intrinsic adducts to neutral loss products. Defaults to False.

    Returns:
        ms2_fragment_dict -- dictionary of sequences and corresponding MS2
                        fragments in format: {
                            seq: {
                                "frag": [masses],
                                "signatures": {
                                    "signature_frag": [masses]
                                    }
                                }
                            }
        where seq = sequence, "frag" = fragment id, masses = MS2 fragment m/z
        values, "signature_frag" = signature fragment id.
    """
    print(f'building MS2 fragment series for {len(sequences)} sequences')
    

    ms2_fragment_dict = {
        sequence: build_multiple_fragment_series_single_sequence(
            sequence=sequence,
            fragment_series=fragment_series,
            adducts=adducts,
            mode=mode,
            add_signatures=add_signatures,
            signatures=signatures,
            losses=losses,
            max_total_losses=max_total_losses,
            loss_product_adducts=loss_product_adducts
        )
        for sequence in sequences
    }

    if uniques:
        print(f'identifying unique fragments for {len(ms2_fragment_dict)} sequences')
        ms2_fragment_dict = generate_unique_fragment_ms2_massdict(
            ms2_fragment_dict)

    return ms2_fragment_dict


def find_unique_fragments_isobaric_set(isobaric_sequencedict):
    """
    Finds unique MS2 fragments from an MS2 fragment mass dictionary in which
    all sequences within the dictionary are isobaric (i.e. have the same
    composition).

    Args:
        isobaric_sequencedict (dict): dictionary of isobaric sequences and
                                    MS2 fragments.

    Returns:
        isobaric_sequencedict: input dict, with list of unique fragments added
                                    to MS2 fragments.
    """
    # initiate list to store fragment masses, for later use in counting
    # occurence of each fragment mass
    all_fragment_masses = []

    # iterate through each isobaric sequence in the input isobaric_sequencedict
    for sequence in isobaric_sequencedict:

        # retrieve fragment dict for sequence
        frag_dict = isobaric_sequencedict[sequence]

        # for each fragment, add associated m/z values to all_fragment_masses
        for fragment, masses in frag_dict.items():

            # do not include signature fragments, as these are a very different
            # beast
            if fragment.find('signature') == -1:
                all_fragment_masses.extend(masses)

    unique_fragment_dict = {}

    # iterate through isobaric sequences again
    for sequence in isobaric_sequencedict:

        # initiate list to store unique fragments, then retrieve fragment dict
        uniques = []
        frag_dict = isobaric_sequencedict[sequence]

        # for each fragment and associated m/z values, count the total occurence
        # of m/z values in the full fragment mass list (all_fragment_masses)
        for fragment, masses in frag_dict.items():
            if fragment.find('signature') == -1:
                mass_count = sum(
                    [all_fragment_masses.count(mass)
                    for mass in masses]
                )

                # check if each fragment m/z value is found only once in the
                # full list of all_fragment_masses; if this is the case, the
                # fragment is unique => add fragment to list of uniques
                if mass_count == len(masses):
                    uniques.append(fragment)

        unique_fragment_dict[sequence] = uniques

    print(f'unique_fragment_dict ={unique_fragment_dict}')
    return unique_fragment_dict


def generate_unique_fragment_ms2_massdict(massdict):
    """
    Takes an MS2 fragment mass dict and identifies unique fragments for each
    sequence in the dict, adds them as a key-value pair to output dict.

    Args:
        massdict (dict): MS2 mass dict in the format: {
            sequence: {
                'frag': [masses]
                }
            }
        where frag = fragment id, masses = m/z values of fragments.

    Returns:
        unique_fragment_dict: same as input massdict, but with extra key-value
                pair added to each sequence subdictionary in format:
                {
                    sequence: {
                        'frag': [masses],
                        'unique_fragments': ['frag_1', 'frag_2'...]
                    }
                }
                where 'frag_1', 'frag_2' = fragment string ids for fragments
                that are unique to sequence.
    """

    # generate dictionary of isobaric sets in format:
    #       {
    #           sorted(seq): [seqs]
    #   }
    # where sorted(seq) = alphabetically sorted sequence; seqs = list of
    # sequences that have same composition (i.e. sorted(seq))
    isobaric_dict = generate_dict_isobaric_sequences(massdict.keys())

    # initiate dict to store ms2_mass_fragment dicts with unique fragments
    # added
    unique_fragment_dict = {}

    # iterate through isobaric sets, updating sequence identifying fragments
    # unique to sequences within each set
    for sequences in isobaric_dict.values():
        sequence_dict = {
            sequence: massdict[sequence] for sequence in sequences
        }

        # update unique_fragment_dict with full MS2 sequence dict for sequence
        # in isobaric_group, with unique fragments added for each sequence
        # within the group
        unique_fragment_dict.update(find_unique_fragments_isobaric_set(sequence_dict))

    # return full MS2 sequence fragment dict, with unique fragments included
    return unique_fragment_dict

def find_fragments_by_index(
    fragments,
    target_index
):
    """
    Takes a list of fragment ids and returns fragments within list whose
    position in the fragment series matches the target_index specified.
    Fragments must be denoted by ONE LETTER codes and index.

    Args:
        fragments (list): list of fragment ids (strings).
        target_index (int): target fragment index.

    Returns:
        fragments (list): list of fragments that match target index (e.g.
                    ['y14', 'b14'] if target_index = 14).
    """
    # initialise fragments list
    fragments_to_mod = []
    
    for fragment in fragments:
        
        # remove one letter code from fragments, leaving only string of fragment
        # index - e.g. 'b14' => '14', 'y2' => '2'   
        fragment_no = int(fragment[1::])
        
        # return list of fragments with index matching target
        if fragment_no == target_index:
            fragments_to_mod.append(fragment)
            
    return fragments_to_mod

def add_terminal_modification_sequence_all_fragments(
    sequence,
    fragment_dict,
    terminal_modifications_dict,
    universal_shift=False
):
    """ This function takes a sequence, all associated subsets of MS2 fragments,
    and adds a terminal modification onto MS2 fragment masses.
    
    Arguments:
        sequence {str} -- full length sequence string comprised of monomer one
                        letter codes
        fragment_dict {dict} -- dictionary of MS2 fragments (different fragment
                        types are allowed as long as they have the same
                        terminus and start, end positions) and associated
                        unmodified masses
        modification_mass {float} -- full neutral monoisotopic mass of modification
        modification_terminus {float} -- mass lost upon addition of modification
                        to fragment mass(es)
    
    Keyword Arguments:
        universal_shift {bool} -- specifies whether every fragment mass
                        must be shifted by modification mass; if True, only
                        modified masses are returned; if False, masses of
                        modified + unmodified fragments are returned.
                        (default: {False})
    
    Returns:
        dict -- format: keys = fragments, values = list of modified and unmodified masses
    """
    # get a list of fragment types
    fragment_series = list(
        set(
            [frag[0] for frag in fragment_dict if frag != 'signature']
        )
    )
    
    # initiate dictionary of ordered fragments
    ordered_fragments = {}
    for frag in fragment_series:
        ordered_fragments[frag] = {}
    
    for fragment, mass in fragment_dict.items():
        
        # initalise fragment, mass dictionary
        frag_mass = {}
        frag_mass[fragment] = mass
        frag_type = fragment[0]
        
        # order fragments by fragment type
        if frag_type in ordered_fragments.keys():
            ordered_fragments[frag_type].update(frag_mass)
    
    # initiate dictionary of all modified fragments 
    modified_fragdict = {}
    
    # for each fragment type generate modified fragments
    for fragment_type, fragment_dict in ordered_fragments.items():
        if fragment_type in fragment_series:
    
            # for each fragment type add terminal modifications to each fragment
    
            # get information for fragment type
            frag_info = FRAG_SERIES[fragment_type]
            terminus = frag_info['terminus']

            # generate modified fragments for each fragment type
            modified_fragments = add_terminal_ms2_modification_sequence_fragment_subsets(
                                 sequence=sequence,
                                 fragment_dict=ordered_fragments[fragment_type],
                                 frag_terminus=terminus,
                                 terminal_modifications_dict=terminal_modifications_dict,
                                 universal_shift=universal_shift)
        
        # add modified fragments to dictionary
        modified_fragdict.update(modified_fragments)
        
    return modified_fragdict

def add_terminal_ms2_modification_sequence_fragment_subsets(
    sequence,
    fragment_dict,
    frag_terminus,
    terminal_modifications_dict,
    universal_shift=False
):
    """
    Takes a sequence, associated subset of MS2 fragments, and adds a terminal
    modification onto MS2 fragment masses. IMPORTANT: ALL FRAGMENT TYPES IN
    fragment_dict MUST HAVE THE SAME TERMINUS, START AND END POSITION.

    Args:
        sequence (str): full length sequence string comprised of monomer one
                        letter codes.
        fragment_dict (dict): dictionary of MS2 fragments (different fragment
                        types are allowed as long as they have the same
                        terminus and start, end positions) and associated
                        unmodified masses.
        frag_terminus (int): 0 or -1 for start and end terminus, respectively;
                        retrieved from FRAG_SERIES['terminus'] in polymer
                        config file.
        modification_mass (float): full neutral monoisotopic mass of modification.
        modification_massdiff (float): mass lost upon addition of modification
                        to fragment mass(es).
        modification_terminus (int): either 0 or -1 for start and end terminus,
                        respectively; specifies whether terminal modification
                        is to be added to start or end terminus.
        universal_shift (bool, optional): specifies whether every fragment mass
                        must be shifted by modification mass; if True, only
                        modified masses are returned; if False, masses of
                        modified + unmodified fragments are returned.
                        Defaults to False.

    Returns:
        modified_fragdict: dictionary of fragment_ids and associated masses,
                        with modification masses added to fragments containing
                        target terminus.
    """
    # initiate dictionary to store fragments with modified masses
    modified_fragdict = {}
    
    for terminal_modification in terminal_modifications_dict.values():
        if terminal_modification:
            # get mass, mass diff and terminus of modification
            modification_mass=MODIFICATIONS[terminal_modification[0]]["mass"]
            modification_massdiff=MODIFICATIONS[terminal_modification[0]]["mass_diff"]["ms1"]
            modification_terminus=MODIFICATIONS[terminal_modification[0]]["termini"]

    # if fragment series begins at the same terminus as the
    # modification_terminus, the terminal modification will be added to all
    # fragments
        if frag_terminus == modification_terminus:
            # iterate through fragment and add terminal modification to
            # fragment masses
            for fragment, masses in fragment_dict.items():
                modified_fragdict[fragment] = add_modification_sequence_mass_list(
                    mass_list=masses,
                    modification_mass=modification_mass,
                    mod_mass_diff=modification_massdiff,
                    universal_shift=universal_shift
                )
            
        # if fragment series begins at opposite terminus to modification_terminus,
        # only fragments that contain modification_terminus (i.e. full-length
        # fragments that span the whole sequence) will be modified
        elif frag_terminus != modification_terminus:
            
            # find fragments that contain target_terminus (i.e. final fragment
            # in series)
            modified_fragments = find_fragments_by_index(
                fragments=fragment_dict.keys(),
                target_index=len(sequence)
            )
            
            # add modification only to fragments that contain the target terminus
            # for modification
            for fragment in modified_fragments:
                modified_fragdict[fragment] = add_modification_sequence_mass_list(
                    mass_list=fragment_dict[fragment],
                    modification_mass=modification_mass,
                    mod_mass_diff=modification_massdiff,
                    universal_shift=universal_shift
                )
                
                # update modified_fragdict with unmodified fragments so that all
                # fragments and associated masses are returned
                modified_fragdict.update(
                    {
                        fragment : masses
                        for fragment, masses in fragment_dict.items()
                        if fragment not in modified_fragments
                }
            )

    return modified_fragdict