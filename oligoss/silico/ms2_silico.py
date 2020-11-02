"""
This file contains functions for generating MS2 silico libraries
"""
from .helpers.helpers import (
    reverse_sequence,
    return_monomer_ids_sequence,
    get_full_neutral_masses_sequence
)

from functools import lru_cache


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

        build_single_fragment.cache_clear()

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
        cache_clear (bool, optional): specifies whether to clear cached data
            for memoized build_single_fragment function
    Returns:
        Dict[str, Dict[str, List[float]]]: dict of sequences and corresponding
            fragment sub dictionaries for fragment series.
    """

    #  build fragments for sequences
    ms2_dict = {
        sequence: build_single_fragment_series_sequence(
            sequence=sequence,
            params=params,
            polymer=polymer,
            series_id=fragment_series
        )
        for sequence in sequences
    }

    return ms2_dict

def build_single_fragment_series_sequence(
    sequence,
    params,
    polymer,
    series_id
):
    """
    Takes a sequence and builds a single MS2 fragment series subdict.

    Args:
        sequence (str): sequence string
        params (Parameters): Parameters object
        polymer (Polymer): Polymer object
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

    #  get range of indices for possible exceptions to apply
    exception_range = find_exception_range(
        sequence_monomer_ids=tuple(sorted(monomer_ids)),
        series_id=series_id,
        polymer=polymer
    )

    #  iterate through sequence backbone, creating subsequences and building
    #  fragments
    for i in range(polymer.fragment_info[series_id]["start"], end):

        # create tuple but leave modification in place if terminal (0)
        if '[' in monomer_ids[0]:
            sorted_subsequence = (
                monomer_ids[0],) + tuple(sorted(monomer_ids[1:i + 1]))
        else:
            sorted_subsequence = tuple(sorted(monomer_ids[0:i + 1]))

        #  call memoized function to get neutral_masses for subsequence
        #  NOTE: for memoization to work, make sure subsequence monomers list
        #  is ALWAYS sorted before being passed into memoized neutral mass func
        neutral_masses = get_memoized_ms2_neutral_masses(
            subsequence_monomer_list=sorted_subsequence,
            params=params,
            polymer=polymer
        )

        #  create subsequence from subsequence monomer id strings
        subsequence_monomers = tuple(
            monomer_ids[polymer.fragment_info[series_id]["start"]:i + 1])

        #  call another memoized function to get m/z values for subsequence
        fragment_mz_list = build_single_fragment(
            exception_range=exception_range,
            subsequence_monomers=subsequence_monomers,
            neutral_masses=tuple(neutral_masses),
            polymer=polymer,
            params=params,
            series_id=series_id
        )

        #  add fragment m/z values to fragment_subdict for sequence
        fragment_subdict[f'{series_id}{i + 1}'] = fragment_mz_list

    return fragment_subdict

@lru_cache(maxsize=10)
def find_exception_range(
    sequence_monomer_ids,
    series_id,
    polymer
):
    """
    Takes a sequence monomer list and works out what range of indices to look
    for possible exceptions to standard MS2 fragmentation rules.

    Args:
        sequence_monomer_ids (Tuple[str]): monomer id strings in sequence /
            composition. NOTE: as this is a memoized function, it is better to
            pass in monomer ids for composition, not unique sequences - this
            will assume isomeric sequences have same cache entry.
        series_id (str): fragment series id.
        polymer (Polymer): polymer object.

    Returns:
        Tuple[int]: (e_start, e_end) where e_start, e_end = start, end index
            in sequence for looking for potential fragmentation exceptions.
    """

    #  there are no polymer exceptions relevant to this fragment series, return
    #  NoneType
    if not polymer.exceptions or series_id not in polymer.exceptions:
        return None

    #  get list of monomer ids for monomers that potentially cause exceptions
    #  for this fragment series. Return NoneType if none are present
    exception_monomers = [
        x for x in sequence_monomer_ids
        if x in polymer.exceptions[series_id]]
    if not exception_monomers:
        return None

    #  keep track of where exceptions start and end
    e_start, e_end = 0, len(sequence_monomer_ids)

    #  iterate through monomers that cause exceptions, keeping track of where
    #  these potential exceptions start and end. Update absolute indices of
    #  e_start and e_end accordingly
    for mon in exception_monomers:
        for j, info in enumerate(polymer.exceptions[series_id][mon].values()):
            if j == 0:
                if info["start"] > 0:
                    e_start = info["start"]
                    e_end = e_end - info["end"]
            else:
                e_start = min(e_start, info["start"])
                e_end = max(e_end, len(sequence_monomer_ids) - info["end"])

    #  finally, return range of indices to look for exceptions in sequence
    return (e_start, e_end)

@lru_cache(maxsize=10000)
def build_single_fragment(
    exception_range,
    subsequence_monomers,
    neutral_masses,
    polymer,
    params,
    series_id
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
        exception_range (Tuple[int]): range of indices at which exceptions to
            standard fragmentation rules are relevant. Format:
                (e_start, e_end) where e_start, e_end = initial and final
                absolute index in subsequence string at which exceptions could
                apply.
        subsequence_monomers (Tuple[str]): str monomer ids for subsequence.
            NOTE: these MUST be in the order they appear in the subsequence.
        neutral_masses (List[float]): list of neutral masses for subsequence
        polymer (Polymer): Polymer object
        params (Parameters): Parameters object
        series_id (str): fragment series id string. NOTE: this is NOT the id
            of the individual fragment (e.g to build a 'b1' or 'b2' fragment,
            series_id == 'b')

    Returns:
        List[float]: list of fragment m/z values.
    """

    #  get mass_diff for fragment series
    mass_diff = polymer.fragment_info[series_id]["mass_diff"]

    #  get intrinsic exchangeable ions in fragment series
    i_adduct = polymer.fragment_info[series_id]["intrinsic_adduct"]

    #  get intrinsic charge for non-exchangeable ions in fragment series
    i_charge = polymer.fragment_info[series_id]["intrinsic_charge"]

    #  intrinsic charge exceeds charge range, can't build fragments
    if i_charge > params.silico.ms2.max_z:
        return []

    #  check whether exceptions to standard fragmentation rules could apply
    if exception_range:
        exceptions = apply_fragmentation_exceptions(
            subsequence_monomers=subsequence_monomers,
            exception_range=exception_range,
            series_id=series_id,
            polymer=polymer
        )
        mass_diff = exceptions["mass_diff"]
        i_adduct = exceptions["intrinsic_adduct"]
        i_charge = exceptions["intrinsic_charge"]

    #  apply fragment mass_diff to neutral masses to get default fragment masses
    default_masses = [mass + mass_diff for mass in neutral_masses]

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

@lru_cache(maxsize=1000)
def add_signature_ions_sequence(
    sequence,
    params,
    polymer
):
    """
    Takes a sequence and returns a dict of monomer-specific signature ions

    Args:
        sequence (str): sequence str
        params (Parameters): Parameters object
        polymer (Polymer): Polymer object

    Returns:
        Dict[str, List[float]]: dict of signature ions and associated m/z values
    """
    #  get signature types for sequence
    signature_types = params.silico.ms2.signatures

    #  get monomer strings from sequence
    monomers = return_monomer_ids_sequence(
        sequence=sequence,
        return_modified=False,
        return_set=True
    )

    #  init dict to store signature ions
    ms2_signatures = {}

    #  iterate through signature types, and if any monomers have signature ions
    #  of this type add them to ms2_signatures dict
    for sig_type in signature_types:
        signatures = polymer.signatures[sig_type]
        for monomer in monomers:
            if monomer in [x[0] for x in signatures]:
                ms2_signatures[f'{sig_type}{monomer}'] = [
                    x[1] for x in signatures
                    if x[0] == monomer
                ]

    return ms2_signatures

def apply_fragmentation_exceptions(
    subsequence_monomers,
    exception_range,
    series_id,
    polymer
):
    """
    Takes a subsequence with its possible fragmentation exceptions and works
    out which exceptions to apply to the following MS2 fragmentation
    properties:
        1. mass_diff
        2. intrinsic_adduct
        3. intrinsic_charge

    Args:
        subsequence_monomers (List[str]): list of monomer strings for
            subsequence. NOTE: these MUST be in the order they appear in the
            subsequence.
        exception_range (Tuple[int]): (e_start, e_end). e_start, e_end =
            absolute index at which exceptions start and end, respectively.
        series_id (str): fragnent series one letter code.
        polymer (Polymer): polymer object.

    Returns:
        dict: dict of fragmentation properties in format:
            {
                "mass_diff": mass_diff (float),
                "intrinsic_charge": intrinsic_charge (float),
                "intrinsic_adduct": intrinsic_adduct (tuple)}
    """
    #  init dict to store values for fragmentation properties (m_diff, i_charge
    #  and i_adduct). For each of these, if exceptions do not apply the default
    #  values will be returned.
    frag_prop_values, frag_index = {}, len(subsequence_monomers)

    fragment_info = polymer.fragment_info[series_id]
    if frag_index <= exception_range[0] or frag_index >= exception_range[1]:
        return fragment_info

    #  list of properties with potential exceptions to standard fragmentation
    exception_props = ["mass_diff", "intrinsic_charge", "intrinsic_adduct"]

    #  iterate through subsequence monomers, check if monomer can cause
    #  fragmentation exceptions and apply as necessary
    for monomer in subsequence_monomers:
        if monomer in polymer.exceptions[series_id]:
            for prop in exception_props:
                if prop in polymer.exceptions[series_id][monomer]:
                    exception_info = polymer.exceptions[
                        series_id][monomer][prop]
                    for position in exception_info["positions"]:
                        if subsequence_monomers[position] == monomer:
                            if prop == "mass_diff":
                                frag_prop_values[prop] = exception_info[
                                    "exception_value"]
                            else:
                                frag_prop_values[prop] = exception_info[
                                    "exception_value"][prop]
                else:
                    frag_prop_values[prop] = fragment_info[prop]
                if prop not in frag_prop_values:
                    frag_prop_values[prop] = fragment_info[prop]

    if not frag_prop_values:
        return fragment_info

    return frag_prop_values
