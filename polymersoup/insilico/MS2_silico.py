"""
This file contains functions for generating theoretical MS2 data for ALife
Polymer Data
"""
from .Constants.GlobalChemicalConstants import *
from .Config_files.Depsipeptide_config import *
from .helpers.helpers import *

def __init__():
    """
    [init]
    """
def build_fragment_series_single_sequence(
    sequence,
    fragment_series,
    adducts,
    mode='pos',
    losses=True,
    max_total_losses=None,
    loss_product_adducts=None
):
    """
    This function takes a sequence and builds an MS2 fragment series,
    returning a dictionary of fragment ids and corresponding m/z values

    Arguments:
        sequence {str} -- sequence string of constituent monomer one letter
                            codes
        fragment_series {str} -- fragment string used to denote fragment
                                series

    Keyword Arguments:
        mode {str} -- specifies whether fragmented species are cationic or
                    anionic, use 'pos' for cations and 'neg' for anions
                    (default: {'pos'})

    Returns:
        fragment_dict -- dictionary of fragments and corresponding m/z
                            values
    """
    # get fragment terminus (i.e. end of sequence fragment series starts from)
    terminus = FRAG_SERIES[fragment_series]['terminus']

    # get fragment mass diff (i.e. mass difference between fragments and
    # equivalent sub-sequences)
    mass_diff = FRAG_SERIES[fragment_series]['mass_diff'][mode]
    if not mass_diff:
        print(f'{fragment_series} fragment series does not exist for {mode} mode')
        return {}

    # get fragment series starting position along linear polymer (i.e. number
    # of residues from terminus to start building fragment series)
    start = FRAG_SERIES[fragment_series]['start']

    # get fragment series ending position along linear polymer (i.e. number of
    # residues from terminus to stop building fragment series)
    end = FRAG_SERIES[fragment_series]['end']

    # fragment series begins at last residue, reverse the sequence before
    # building the series (e.g. peptide 'y' fragments)
    if terminus == -1:
        sequence = sequence[::-1]

    # initiate fragment dictionary to be populated by fragments and associated
    # masses
    fragment_dict = {}

    # build fragment dict in format { f'{frag_series}{i}' : mass} where
    # frag_series = fragment series code (str - one character for
    # non-signatures)

    # check whether to include side chain-specific loss products in fragment
    # mass lists; if not, fragment_dict is returned with full fragment masses
    # only

    # initiate dict to store fragment side chain-specific loss product masses
    loss_fragment_dict = {}

    # iterate through fragments, building and adding to fragment dict one at a
    # time
    for i in range(start, len(sequence)-end):
        fragment, sub_sequence = f'{fragment_series}{i+1}', sequence[0:i+1]
        fragment_mass = find_sequence_mass(sub_sequence) + mass_diff
        fragment_dict[fragment] = [float(f'{fragment_mass:.4f}')]
        if losses:
            loss_fragment_dict[fragment] = add_sidechain_neutral_loss_products_sequence(
                sub_sequence,
                fragment_mass,
                max_total_losses
            )

    # add adducts to full fragment masses (non-loss products)
    fragment_dict = add_adducts_ms2_fragments(
        fragment_dict,
        adducts,
        mode,
        fragment_series)

    # add adducts to ms2 loss products if any are specified for loss products
    if loss_product_adducts:
        loss_fragment_dict = add_adducts_ms2_fragments(
            loss_fragment_dict,
            loss_product_adducts,
            mode,
            fragment_series)

    # combine fragment dict with loss product fragment dict
    for fragment, masses in loss_fragment_dict.items():
        masses.extend(fragment_dict[fragment])
        masses = sorted(list(set(masses)))
        fragment_dict[fragment] = masses

    # finally, return MS2 fragment dictionary
    return fragment_dict

def build_multiple_fragment_series_single_sequence(
    sequence,
    fragment_series,
    adducts,
    mode='pos',
    add_signatures=True,
    signatures=None,
    losses=True,
    max_total_losses=None,
    loss_product_adducts=None
):
    """
    Takes a single sequence, list of fragment series (e.g. ['b', 'y'] for
    standard peptide b- and y- fragments) and outputs a dictionary of MS2
    fragments and corresponding m/z values

    Args:
        sequence (str): Sequence string comprised of monomer one letter codes
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
        fragment_dict: dictionary of MS2 fragments and corresponding m/z values
    """
    # if no fragment series are specified, all fragment types in FRAG_SERIES
    # dict from polymer config file are to be incorporated into MS2 fragment
    # dict
    if not fragment_series:
        fragment_series = FRAG_SERIES.keys()

    # initiate fragment dictionary to be populated
    fragment_dict = {}

    # iterate through fragment series list and update fragment dictionary for
    # each fragment type
    for series in fragment_series:
        fragment_dict.update(
            build_fragment_series_single_sequence(
                sequence,
                series,
                adducts,
                mode,
                losses,
                max_total_losses,
                loss_product_adducts))

    # checks if monomer-specific MS2 signature ion fragments are to be added
    # if specified, signature ions are then added to fragment dict
    if add_signatures:
        fragment_dict.update({
            'signatures': generate_monomer_ms2_signature_ions_sequence(
                sequence,
                signatures)})
    return fragment_dict

def add_adducts_ms2_fragments(
    fragment_dict,
    adducts,
    mode='pos',
    fragment_series=None,
    return_base_fragments=True
):
    """
    This function takes a dictionary of MS2 fragments and associated m/z
    values, and adds singly charged adducts

    Arguments:
        fragment_dict {dict} -- dictionary of fragments and associated m/z
                            values. keys = fragment id strings, values = lists
                            of associated masses. See function:
                            build_fragment_series_single_sequence
        adducts {list} -- list of adducts to be added to MS2 fragments

    Keyword Arguments:
        mode {str} -- either 'pos' or 'neg' for positive and negative mode,
                        respectively (default: {'pos'})
        fragment_series {list} -- list one character strings denoting fragment
                        series. If None, adducts will be added to all fragment
                        series that do not correspond to signature ions
                        (default: {None})
        return_base_fragments {bool} -- specific whether to return input
                        fragment masses in adduct list (default: {True})

    Raises:
        Exception: raised if any of the input adducts do not match the
                    detection mode (i.e. anions in positive mode, cations in
                    positive mode). Functions for complex adducts are yet to
                    be written

    Returns:
        adduct_fragdict -- dictionary of fragment string codes and associated
                        m/z value lists with adducts. ALL IN +1 CHARGE STATE
    """
    # list fragments, excluding monomer-specific signature ion fragments
    fragments = [key for key in fragment_dict if key != 'signatures']

    # get list of fragments that belong to specific fragment series,
    # if specified
    if fragment_series:
        fragments = [
            fragment[0] for fragment in fragments
            if fragment[0] in fragment_series]

    # list unique types of fragments - this is denoted by the first character
    # of fragment id strings (i.e. fragment one letter codes)
    fragment_series = list(set([fragment[0] for fragment in fragments]))

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
                print(f'type mass, i_adduct = {type(mass)}, {type(intrinsic_adduct)}')
                mass -= intrinsic_adduct
                adduct_masses.extend(add_adducts_sequence_mass(
                    mass,
                    adducts,
                    1,
                    1,
                    mode))

            # if standard fragments (i.e. with unchanged input masses) are to
            # be returned, include these in output adduct_masses list
            if return_base_fragments:
                adduct_masses.extend(masses)

            # add fragment and associated adducts to output adduct_fragdict,
            # removing any duplicate masses - which should not be there anyway
            fragment_adducts[fragment] = sorted(list(set(adduct_masses)))

    # finally, return dictionary of fragments and associated adduct mass lists
    return fragment_adducts

def generate_monomer_ms2_signature_ions_sequence(
    sequence,
    signatures=None,
):
    """
    This function generates a dictionary of monomer signature ions for
    a sequence
    Arguments:
        sequence {str} -- sequence string with constituent monomer one letter
                            codes
    Keyword Arguements:
        signatures {list} -- list of signature fragment codes - e.g. ['Im']
                            for amino acid immonium fragments. If None, all
                            signature ion types will be returned in signature
                            ion dictionary
                            (default: {None})
    Returns:
        MS2_signature_dict -- dictionary of monomer signature fragments
                            and associated m/z values
    """
    # get list of unique monomers within sequence
    monomers = list(set([c for i, c in enumerate(sequence)]))

    # get list of signature ion id strings used to access signature ion values
    # from polymer config file
    if not signatures:
        signatures = MS2_SIGNATURE_IONS.keys()

    # initiate output MS2 signature ion dict
    MS2_signature_dict = {}

    # iterate through signature ion types and add to signature ion dictionary
    for signature in signatures:
        signature_ions = MS2_SIGNATURE_IONS[signature]

        MS2_signature_dict.update(
            {
                f'{signature}{monomer}': signature_ions[monomer]
                for monomer in monomers
                if monomer in signature_ions
            }
        )

    # finally, return dictionary of monomer signature ions and associated
    # lists of m/z values
    return MS2_signature_dict

def generate_ms2_mass_dictionary(
    sequences,
    fragment_series,
    adducts,
    mode,
    add_signatures=True,
    signatures=None,
    losses=True,
    max_total_losses=None,
    loss_product_adducts=None
):
    """
    Takes a list of sequences and generates full ms2 fragment dictionary

    Args:
        sequences (str): sequence string comprised of monomer one letter codes
        fragment_series (list): list of fragment types - must be keys in
                            FRAG_SERIES dict in polymer config file
        adducts (list): list of adduct strings
        mode (str): either 'pos' or 'neg' for positive and negative mode,
                    respectively
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
        loss_product_adducts (list, optional): list of adducts to add to loss
                            product fragment masses. Defaults to None.

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
        values, "signature_frag" = signature fragment id
    """
    ms2_fragment_dict = {
        sequence: build_multiple_fragment_series_single_sequence(
            sequence,
            fragment_series,
            adducts,
            mode,
            add_signatures,
            signatures,
            losses,
            max_total_losses)
        for sequence in sequences}

    return ms2_fragment_dict
