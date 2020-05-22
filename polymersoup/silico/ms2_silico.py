"""
This file contains functions for generating MS2 silico libraries
"""
from functools import lru_cache

from ..utils.global_chemical_constants import (
    ANIONS,
    CATIONS
)
from .helpers.helpers import (
    reverse_sequence,
    return_monomer_ids_sequence,
    get_full_neutral_masses_sequence
)

def build_linear_fragments_sequence_dict(
    sequences,
    params,
    polymer
):
    """
    Builds full linear MS2 fragment dict for sequences

    Args:
        sequences (List[str]): list of sequence strings
        params (Parameters): Parameters object
        polymer (Polymer): Polymer object

    Returns:
        Dict[str, Dict[str, List[float]]]: Dict of sequences and corresponding
            linear MS2 fragment subdicts in format: {
                seq: {
                    frag_id: [m/z],
                    ...
                },
                ...
            }
            where seq = sequence, frag_id = fragment id, [m/z] = list of m/z
            values for individual MS2 fragments
    """

    #  init dict to store sequences, fragment subdicts
    ms2_dict = {}

    #  iterate through linear fragment series, build fragment dicts for each
    for series in polymer.fragment_info:
        series_dict = build_fragment_dict_for_sequences(
            sequences=sequences,
            params=params,
            polymer=polymer,
            fragment_series=series
        )
        if not ms2_dict:
            ms2_dict = series_dict
        else:
            for sequence, frag_dict in ms2_dict.items():
                frag_dict.update(series_dict[sequence])

    #  clear cache for memoized function to save memory
    get_memoized_ms2_neutral_masses.cache_clear()

    return ms2_dict

def build_fragment_dict_for_sequences(
    sequences,
    params,
    polymer,
    fragment_series
):
    """
    Takes a list of sequences and builds MS2 fragment dict for a single fragment
    series.

    Args:
        sequences (List[str]): list of sequence strings
        params (Parameters): Parameters object
        polymer (Polymer): Polymer object
        fragment_series (str): string id for a single fragment series

    Returns:
        Dict[str, Dict[str, List[float]]]: dict of sequences and corresponding
            fragment sub dictionaries for fragment series.
    """
    #  get fragment intrinsic adduct
    i_adduct = retrieve_intrinsic_adduct_series(
        fragment_series=fragment_series,
        polymer=polymer,
        params=params
    )

    #  build fragments for sequences
    ms2_dict = {
        sequence: build_single_fragment_series_sequence(
            sequence=sequence,
            params=params,
            polymer=polymer,
            i_adduct=i_adduct,
            series_id=fragment_series
        )
        for sequence in sequences
    }

    #  clear cache for memoized fragment builder function as next set of calls
    #  won't apply to this fragment series
    build_single_fragment.cache_clear()

    return ms2_dict

def build_single_fragment_series_sequence(
    sequence,
    params,
    polymer,
    i_adduct,
    series_id
):
    """
    Takes a sequence and builds a single MS2 fragment series subdict.

    Args:
        sequence (str): sequence string
        params (Parameters): Parameters object
        polymer (Polymer): Polymer object
        i_adduct (List[float]): [mass, min_z, max_z] for fragment series'
            intrinsic adduct
        series_id (str): fragment series id string

    Returns:
        Dict[str, List[float]]: dict of individual fragment ids and
            corresponding list of m/z values.
    """

    #  init dict to store fragments and m/z values for sequence / series
    fragment_subdict = {}

    #  reverse sequence if fragment series is indexed from -1 terminus
    if polymer.fragment_info[series_id]["terminus"] == -1:
        sequence = reverse_sequence(sequence=sequence)

    #  retrieve monomer_ids, in order they appear in sequence / reversed
    #  sequence
    monomer_ids = return_monomer_ids_sequence(
        sequence=sequence,
        return_modified=True,
        return_set=False
    )

    #  get absolute ending index for fragment series
    sequence_length = len(monomer_ids)
    end = sequence_length - polymer.fragment_info[series_id]["end"]

    #  iterate through sequence backbone, creating subsequences and building
    #  fragments
    for i in range(polymer.fragment_info[series_id]["start"], end):

        #  sort subsequence monomers so isomeric subsequences have same cache
        #  entry - see below. Convert to make hashable for memoization
        sorted_subsequence = tuple(sorted(monomer_ids[0:i + 1]))

        #  call memoized function to get neutral_masses for subsequence
        #  NOTE: for memoization to work, make sure subsequence monomers list
        #  is ALWAYS sorted before being passed into memoized neutral mass func
        neutral_masses = get_memoized_ms2_neutral_masses(
            subsequence_monomer_list=sorted_subsequence,
            params=params,
            polymer=polymer
        )

        for neutral_mass in neutral_masses:
            if neutral_mass < 70:
                raise Exception(neutral_mass, sorted_subsequence)
        #  create subsequence from subsequence monomer id strings
        subsequence = "".join(
            monomer_ids[polymer.fragment_info[series_id]["start"]:i + 1])

        #  call another memoized function to get m/z values for subsequence
        fragment_mz_list = build_single_fragment(
            subsequence=subsequence,
            neutral_masses=tuple(neutral_masses),
            polymer=polymer,
            params=params,
            series_id=series_id,
            i_adduct=tuple(i_adduct)
        )

        #  add fragment m/z values to fragment_subdict for sequence
        fragment_subdict[f'{series_id}{i + 1}'] = fragment_mz_list

    return fragment_subdict

@lru_cache(maxsize=100000)
def build_single_fragment(
    subsequence,
    neutral_masses,
    polymer,
    params,
    series_id,
    i_adduct
):
    """
    Builds a single fragment of a single type for a particular single
    subsequence.
    NOTE: This is a memoized function as many sequences have
    overlapping subsequences. For most efficient use, it should be called
    successively to build fragments of the SAME type (i.e. same series_id).
    This will ensure the maximum hits for previous cache entries when building
    MS2 fragment dicts. NOTE: please clear cache after building MS2 sequence
    dict for target fragment series as this will save memory.

    Args:
        subsequence (str): subsequence string
        neutral_masses (List[float]): list of neutral masses for subsequence
        polymer (Polymer): Polymer object
        params (Parameters): Parameters object
        series_id (str): fragment series id string. NOTE: this is NOT the id
            of the individual fragment (e.g to build a 'b1' or 'b2' fragment,
            series_id == 'b')
        i_adduct (List[float]): [mass, min_z, max_z] for intrinsic exchangeable
            ion for fragment series

    Returns:
        List[float]: list of fragment m/z values.
    """
    #  get mass_diff for fragment series
    mass_diff = polymer.fragment_info[series_id]["mass_diff"][params.mode]

    #  get intrinsic charge for non-exchangeable ions in fragment series
    i_charge = 0
    if "intrinsic_charge" in polymer.fragment_info[series_id]:
        if polymer.fragment_info[series_id]["intrinsic_charge"][params.mode]:
            i_charge += polymer.fragment_info[
                series_id]["intrinsic_charge"][params.mode]

    #  intrinsic charge exceeds charge range, can't build fragments
    if i_charge > params.silico.ms2.max_z:
        return []

    #  apply fragment mass_diff to neutral masses to get default fragment masses
    default_masses = [
        mass + mass_diff for mass in neutral_masses
    ]

    #  init list to store m/z values for MS2 fragments
    ms2_ions = []

    #  iterate through ms2 adducts and add to fragments
    for adduct_info in params.silico.ms2.adducts.values():
        if i_adduct and i_adduct == adduct_info:
            pass
        else:
            ms2_ions.extend(generate_adduct_for_ms2_fragment(
                params=params,
                adduct_info=adduct_info,
                default_masses=default_masses,
                i_charge=i_charge,
                i_adduct=i_adduct
            ))

    #  if intrinsic adducts, make sure these are kept and included in
    #  fragment m/z list
    if i_adduct[2] and i_adduct[2] > 1:
        ms2_ions.extend(generate_adduct_for_ms2_fragment(
            params=params,
            adduct_info=i_adduct,
            default_masses=default_masses,
            i_charge=i_charge,
            i_adduct=i_adduct
        ))

    #  singly charged intrinsic adduct by default, so return these ions
    elif i_adduct[1] and i_adduct[1] >= params.silico.ms2.min_z:
        ms2_ions.extend(default_masses)

    #  fragments are intrinsically charged within charge range, so return these
    #  ions
    if i_charge in range(params.silico.ms2.min_z, params.silico.ms2.max_z + 1):
        ms2_ions.extend(default_masses)

    return list(set([round(mz, 4) for mz in ms2_ions]))


def generate_adduct_for_ms2_fragment(
    params,
    adduct_info,
    default_masses,
    i_charge,
    i_adduct
):
    """
    Creates charged adducts from MS2 fragments

    Args:
        params (Parameters): Parameters object
        adduct_info (List[float]): [adduct_mass, min_z, max_z] for incoming
            exhangeable ion to be added to fragment
        default_masses (List[float]): default masses (NOTE: not m/z values) for
            fragment
        i_charge (int): intrinsic charge of fragment due to non-exchangeable
            ions
        i_adduct (List[float]): [adduct_mass, min_z, max_z] for intrinsic
            exchangeable ion for fragment

    Returns:
        List[float], optional: list of m/z values for fragment if fragment
            can be ionized under constraints set by Parameters object.
    """
    if adduct_info[1] > params.silico.ms2.max_z:
        return []

    #  get minimum charge state for final adducts, return nothing if this
    #  exceeds max available charge
    min_charge = max(params.silico.ms2.min_z, adduct_info[1] + i_charge)
    if min_charge > params.silico.ms2.max_z:

        return []

    #  make sure to subtract intrinsic adduct when exchanging for incoming
    #  adduct
    adduct_mass = adduct_info[0]
    if i_adduct[0]:
        adduct_mass -= i_adduct[0]

    #  init list to store m/z values for fragment
    ms2_ions = []

    #  get maximum charge => whatever is lowest between adduct and pre-set
    #  experiment max charge in Parameters object
    max_charge = min(params.silico.ms2.max_z, adduct_info[2] + i_charge)

    #  iterate through available charge states and work out final m/z list
    for i in range(min_charge, max_charge + 1):
        ms2_ions.extend([
            (mass + adduct_mass) / (i + i_charge)
            for mass in default_masses
        ])

    return ms2_ions

@lru_cache(maxsize=10000)
def get_memoized_ms2_neutral_masses(
    subsequence_monomer_list,
    params,
    polymer
):
    """
    Gets neutral MS2 masses for a subsequence. NOTE: this is a memoized
    funcion. For memoization to work effectively when building MS2 fragment
    dicts make sure to always pass in subsequence_monomer_list as a sorted list.
    This is to ensure that isomeric subsequences with same neutral masses have
    the same cache entry!

    Args:
        subsequence_monomer_list (Tuple[str]): Tuple of monomer id strings
        params (Parameters): Parameters object
        polymer (Polymer): Polymer object

    Returns:
        List[float]: list of neutral masses, including neutral loss products
            for subsequence
    """

    #  create sequence string from subsequence monomers
    subsequence = "".join(list(subsequence_monomer_list))

    #  return full neutral masses for subsequence
    return get_full_neutral_masses_sequence(
        sequence=subsequence,
        params=params,
        polymer=polymer,
        ms_level=2
    )

def retrieve_intrinsic_adduct_series(
    fragment_series,
    polymer,
    params
):
    """
    Checks for any intrinsic adducts for fragment series.

    Args:
        fragment_series (str): fragment series string id
        polymer (Polymer): Polymer object
        params (Parameters): Parameters object

    Returns:
        (float, float, float, optional): intrinsic adduct mass, minimum charge,
            maximum charge
    """
    if "intrinsic_adducts" not in polymer.fragment_info[fragment_series]:
        return None, None, None

    if params.mode == "pos":
        ions = CATIONS
    else:
        ions = ANIONS

    if polymer.fragment_info[
            fragment_series]["intrinsic_adducts"][params.mode]:
        ion = polymer.fragment_info[
            fragment_series]["intrinsic_adducts"][params.mode]
        return ions[ion]

    return None, None, None
