import itertools
from functools import lru_cache

from ...utils.global_chemical_constants import FUNCTIONAL_GROUPS
from ...utils.silico_utils.silico_utils import SILICO_MOD_MARKERS

def find_sequence_mass(
    sequence,
    polymer,
    n_rounded=6
):
    """
    Takes a sequence string and returns its neutral monoisotopic mass

    Args:
        sequence (str): sequence string comprised of monomer one-letter codes,
            sidechain and terminal modification substrings
        polymer (Polymer): instance of Polymer class
        n_rounded (int, optional): number of decimal places to round sequence
            mass down to. Defaults to 6.

    Returns:
        float: neutral monoisotopic mass of sequence
    """
    #  retrieve list of monomer ids in sequence, in order in which they apear
    #  in sequence string
    monomer_ids = return_monomer_ids_sequence(
        sequence=sequence,
        return_modified=True,
        return_set=False
    )

    #  sum mass of monomers in sequence
    sequence_mass = sum([
        find_monomer_mass(monomer_id=monomer_id, polymer=polymer)
        for monomer_id in monomer_ids])

    if type(polymer.mass_diff) == str:
        try:
            polymer.mass_diff = float(polymer.mass_diff)
        except ValueError:
            if polymer.mass_diff[0] == "-":
                polymer.mass_diff = -FUNCTIONAL_GROUPS[polymer.mass_diff[1::]]
            else:
                polymer.mass_diff = FUNCTIONAL_GROUPS[polymer.mass_diff]

    #  subtract mass diffs from sequence mass
    sequence_mass -= polymer.mass_diff * (len(monomer_ids) - 1)

    return round(sequence_mass, n_rounded)

@lru_cache(maxsize=20)
def find_monomer_mass(monomer_id, polymer):
    """
    Takes a monomer id string and returns neutral monoisotopic mass of monomer

    Args:
        monomer_id (str): monomer id string, including any sidechain
            modifications
        polymer (Polymer): instance of Polymer class

    Returns:
         float: neutral monoisotopic mass of monomer
    """

    #  get appropriate Monomer object from Polymer object
    monomer = list(filter(
        lambda x: x.id == monomer_id[0], polymer.monomers))[0]

    #  if monomer_id has no extra characters it is unmodified, so return
    #  standard neutral mass for unmodified monomer
    if len(monomer_id) == 1:
        return monomer.neutral_mass

    # get list of monomer modifications from input monomer string
    monomer_mods = remove_substrings_sequence(
        sequence=monomer_id,
        return_substrings=True
    )

    #  get neutral monoisotopic masses of covalent modifications attached to
    #  monomer sidechain
    modification_masses = [
        find_modification_mass(
            modification=mod[1:len(mod) - 1],
            polymer=polymer,
            subtract_massdiff=True
        )
        for mod in monomer_mods
    ]

    #  return sum of modifications and monomer as final mass
    return sum(modification_masses) + monomer.neutral_mass

@lru_cache(maxsize=20)
def find_modification_mass(modification, polymer, subtract_massdiff):
    """
    Takes a modification string and returns its corresponding mass.

    Args:
        modification (str): modification string.
        polymer (Polymer): Polymer object.
        subtract_massdiff (bool): specifies whether to return mass with
            modification massdiff subtracted.

    Returns:
        float: modification neutral monoisotopic mass
    """

    mod_obj = get_mod_info_from_id_string(
        modification_str=modification,
        polymer=polymer
    )

    #  set mass, mass_diff as defined in Modification object
    mass, mass_diff = mod_obj.mass, mod_obj.mass_diff["ms1"]

    #  make sure mass_diff is supplied as numeric value
    if type(mass_diff) == str:
        mass_diff = FUNCTIONAL_GROUPS[mod_obj.mass_diff["ms1"]]

    #  return mass, subtracting massdiff is specified
    if subtract_massdiff:
        return mass - mass_diff

    return mass

def get_mod_info_from_id_string(
    modification_str,
    polymer
):
    """
    Returns Modification object with id matching modification string.

    Args:
        modification_str (str): modification id string
        polymer (Polymer): Polymer object

    Raises:
        Exception: raised if modification_str does not match any Modification
            objects in polymer

    Returns:
        Modification: Modification object
    """

    #  iterate through Modification objects in polymer, return object with id
    #  matching modification string
    for modification_list in polymer.modifications.values():
        mod_obj = list(filter(
            lambda mod: mod.id == modification_str,
            modification_list
        ))
        if mod_obj:
            return mod_obj[0]

    raise Exception(
        f'No modifications with id "{modification_str}" are present\n'
        'in Polymer object'
    )

def return_core_sequence(sequence):
    """
    Takes a sequence string with or without terminal and covalent modifications
    and returns the equivalent unmodified or "core" sequence.

    Args:
        sequence (str): sequence string

    Returns:
        str: unmodified sequence string.
    """
    #  remove terminal modifications from sequence
    sequence = remove_terminal_modifications_sequence_string(sequence)

    #  remove sidechain modifications from sequence and return
    return remove_sidechain_modifications_sequence_string(sequence)


def return_monomer_ids_sequence(
    sequence,
    return_modified,
    return_set
):
    """
    Takes a sequence comprised of monomer one-letter codes + / - side chain
    and terminal modification strings, and returns list of monomer ids in the
    sequence. NOTE: monomer id returned can either be exclusively one letter
    codes or include any sidechain modifications (see args: return_modified).

    Args:
        sequence (str): sequence string comprised of monomer one-letter codes
            and / or sidechain and terminal modifications
        return_modified (bool, optional): specifies whether to return monomers
            with associated sidechain modification strings. Defaults to True.
        return_set (bool, optional): specifies whether to return list of
            unique monomers (i.e. list of a set) or not. Defaults to True.

    Returns:
        list: list of monomer ids within input sequence string.
    """
    # check whether to include sidechain modifications in returned monomer
    # list. if not, return only 'core' unmodified monomers
    if not return_modified:
        core_sequence = return_core_sequence(
            sequence=sequence
        )

        #  remove duplicate monomers if set is to be returned
        if return_set:
            return list(set([c for i, c in enumerate(core_sequence)]))

        #  return standard list of monomer ids, in order they appear in sequence
        #  (including duplicate monomer ids)
        return [c for i, c in enumerate(core_sequence)]

    #  return monomer ids with modification substrings (if any)
    return get_modified_monomer_ids(
        sequence=sequence,
        return_set=return_set)

def get_modified_monomer_ids(
    sequence,
    return_set
):
    """
    Takes a sequence string and returns list of monomer ids (including
    associated sidechain modifications).

    Args:
        sequence (str): sequence string with +/- sidechain and terminal
            modifications
        return_set (bool): specifies whether to return list of unique monomers.
            if False, full list of monomers (including duplicates) is returned
            sorted by the order they appear in sequence string. If True,
            list of set is returned (i.e. duplicates removed).

    Returns:
        List[str]: list of monomer id strings
    """
    #  remove terminal modifications from sequence string, leaving only monomers
    sequence = remove_terminal_modifications_sequence_string(sequence)

    # init list to store monomers
    monomers = []

    # init list to store indices in sequence which are part of modification
    # substrings
    mod_indices = []

    # iterate through sequence, retrieving sidechain-modified monomers and
    # adding them to list of monomers
    for i, c in enumerate(sequence):

        if c in SILICO_MOD_MARKERS:

            mod_end = min([
                j for j, d in enumerate(sequence)
                if d == SILICO_MOD_MARKERS[c] and j > i
            ])

            # append monomer and its index to monomers list in format [m, i]
            # where m = monomer string and i = index of its first character in
            # sequence
            monomers.append([sequence[i - 1:mod_end + 1], i - 1])

            # get indices of all characters containing modified monomer
            # substrings so that unmodified monomers can be identified
            mod_indices.extend([
                x for x in range(i - 1, mod_end + 1)
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

def remove_terminal_modifications_sequence_string(sequence):
    """
    Takes a sequence string with terminal modifications and removes characters
    corresponding to terminal modifications. NOTE: sidechain modifications
    remain.

    Args:
        sequence (str): sequence string comprised of monomer one-letter codes,
            sidechain and terminal modification substrings.

    Returns:
        str: sequence string with terminal modification substrings removed
    """

    # check whether sequence has 0 terminal modification; if so, remove
    if sequence[0] == SILICO_MOD_MARKERS['terminal'][0]:

        mod_end = min([
            i for i, c in enumerate(sequence)
            if c == SILICO_MOD_MARKERS['terminal'][1]
        ])

        sequence = sequence[mod_end + 1::]

    # check whether sequence has -1 terminal modification; if so, remove
    if sequence[-1] == SILICO_MOD_MARKERS['terminal'][-1]:

        mod_start = max(
            i for i, c in enumerate(sequence)
            if c == SILICO_MOD_MARKERS['terminal'][0]
        )

        sequence = sequence[0:mod_start]

    return sequence

def remove_sidechain_modifications_sequence_string(sequence):
    """
    Takes a sequence string with sidechain modifications and removes characters
    corresponding to sidechain modifications.
    Args:
        sequence (str): sequence string comprised of monomer one-letter codes,
            sidechain and sidechain modification substrings.

    Returns:
        str: sequence string with sidechain modification substrings removed
    """

    #  get list of indices for characters that fall within a modification
    #  substring, remove each from string
    mod_indices = []
    for i, c in enumerate(sequence):

        if c in SILICO_MOD_MARKERS:

            mod_end = min(
                [
                    j for j, d in enumerate(sequence)
                    if d == SILICO_MOD_MARKERS[c] and j > i
                ]
            )
            mod_indices.extend([j for j in range(i, mod_end + 1)])

    #  return core sequence
    return "".join(
        [
            c for i, c in enumerate(sequence)
            if i not in mod_indices
        ]
    )


def reverse_sequence(sequence):
    """
    Takes a sequence comrpised of monomer one-letter codes and / or any
    terminal and sidechain modification substrings and outputs the reverse
    of that sequence.
    Example:
        '[Ole]S(Trt)GAVS(Trt)' -> 'S(Trt)VAGS(Trt)[Ole]'

    Args:
        sequence (str): sequence string comprised of monomer one letter codes
            and / or modification substrings

    Returns:
        str: reverse of input sequence
    """
    reversed_sequence = sequence[::-1]

    # reverse NON-TERMINAL mod markers to deal with reversed sequence
    r_mod_markers = {
        v: k
        for k, v in SILICO_MOD_MARKERS.items()
        if k != 'terminal'
    }

    #  reverse terminal markers
    r_terminal_markers = {
        SILICO_MOD_MARKERS['terminal'][1]: SILICO_MOD_MARKERS['terminal'][0]
    }

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

        #  iterate through each reversed side chain modification and ensure
        #  the now-reversed modification substring is replaced with the original
        #  non-reversed string
        for r_mod in reversed_mods:
            reversed_sequence = reversed_sequence.replace(r_mod, r_mod[::-1])

            #  find positons of side chain modification substrings in reversed
            #  sequence along with their associated monomers
            for i, c in enumerate(sequence):
                if c == r_mod[::-1][0]:

                    if sequence[i:i + len(r_mod)] == r_mod[::-1]:
                        r_mod_monomer = sequence[i - 1:i + len(r_mod)]
                        non_r_mod_monomer = (
                            f"{r_mod_monomer[1::]}{r_mod_monomer[0]}")
                        mod_monomers[r_mod_monomer] = non_r_mod_monomer

    #  move position of sidechain modification substrings in reversed sequence
    #  to ensure they remain associated with correct monomers
    if mod_monomers:
        for monomer, non_r_monomer in mod_monomers.items():
            reversed_sequence = reversed_sequence.replace(
                non_r_monomer, monomer)

    # check whether any terminal mods have been found in reversed sequence
    # if so, switch them back to regular non-reversed substrings
    if reversed_terminal_mods:

        for r_t_mod in reversed_terminal_mods:
            reversed_sequence = reversed_sequence.replace(
                r_t_mod, r_t_mod[::-1])

    return reversed_sequence

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
        substring_markers (dict, optional): dictionary of chars marking the
            start of substrings and chars marking their end. Defaults to {
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

            mod_indices.extend([x for x in range(i, mod_end + 1)])
            substring = sequence[i:mod_indices[-1] + 1]
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

@lru_cache(maxsize=10000)
def ionize_sequence_precursors(sequence, params, polymer):
    """
    NOTE: Memoized function. If called multiple times for MS1 composition
    strings rather than unique sequence strings, this function should be
    faster.

    Takes a sequence string, Parameters and Polymer objects, and returns list
    of MS1 precursor m/z values

    Args:
        sequence (sequence string): sequence id string
        params (Parameters): Parameters object
        polymer (Polymer): Polymer object

    Returns:
        List[float]: list of m/z values for sequence precursors (MS1 ions)
    """

    #  get full list of neutral monoisotopic masses for sequence - standard
    #  neutral mass + neutral loss products
    neutral_masses = get_full_neutral_masses_sequence(
        sequence=sequence,
        params=params,
        polymer=polymer,
        ms_level=1
    )

    #  generate ions from neutral masses and return, sorted in descending order
    #  by m/z
    return sorted(add_adducts_sequence_masses(
        sequence=sequence,
        masses=neutral_masses,
        params=params,
        polymer=polymer,
        ms_level=1
    ), reverse=True)


def get_full_neutral_masses_sequence(sequence, params, polymer, ms_level):
    """
    Takes a sequence string, Parameters and Polymer object and returns list of
    neutral monoisotopic masses for sequence - including full neutral mass
    and potential loss products

    Args:
        sequence (str): sequence id string
        params (Parameters): Parameters object
        polymer (Polymer): Polymer object
        ms_level (int): level of MS data sequence is to eventually be screened
            for in - this parameter is specified because different neutral loss
            fragmentation limits are placed on MS1 and MS2 data.

    Returns:
        List[float]: list comprised of full neutral monositopic mass of sequence
            plus neutral loss products
    """
    #  get the neutral monoisotopic mass of sequence
    neutral_masses = [find_sequence_mass(
        sequence=sequence,
        polymer=polymer,
        n_rounded=6
    )]

    #  generate list of neutral losses for sequence, add to neutral masses
    neutral_masses.extend(add_sidechain_neutral_loss_products_sequence(
        sequence=sequence,
        sequence_mass=neutral_masses[0],
        polymer=polymer,
        params=params,
        ms_level=ms_level
    ))

    return neutral_masses

def add_adducts_sequence_masses(
    sequence,
    masses,
    polymer,
    params,
    ms_level
):
    """
    Takes a sequence, associated neutral masses, and adds ions to generate
    final m/z values.

    Args:
        sequence (str): sequence string
        masses (List[float]): list of neutral monoisotopic masses for sequence
        polymer (Polymer): Polymer object
        params (Parameters): Parameters object
        ms_level (int): level of MS data sequence is to eventually be screened
            for in - this parameter is specified because different adducts may
            be specified for MS1 and MS2.

    Raises:
        Exception: raised if params.mode not in ["pos", "neg"].

    Returns:
        List[float]: list of final m/z values for sequence.
    """

    #  get minimum and max charge set by instrument parameters
    if ms_level == 1:
        min_z, max_z = params.silico.ms1.min_z, params.silico.ms1.max_z
    else:
        min_z, max_z = params.silico.ms2.min_z, params.silico.ms2.max_z

    if ms_level == 1:
        adducts = params.silico.ms1.adducts
    else:
        adducts = params.silico.ms2.adducts

    #  init list to store final m/z values, iterate through adducts to extend
    final_mz = []
    for adduct in adducts:
        if min_z:
            min_charge = max(min_z, adducts[adduct][1])
        else:
            min_charge = adducts[adduct][1]
        if max_z:
            max_charge = min(max_z, adducts[adduct][2])
        else:
            max_charge = adducts[adduct][2]
        if min_charge <= max_charge:
            for x in range(min_charge, max_charge + 1):
                final_mz.extend([
                    round((mass + adducts[adduct][0]) / x, 4) for mass in masses
                ])

    return final_mz

@lru_cache(maxsize=100000)
def retrieve_loss_products_sequence_monomers(monomer_ids, polymer):

    #  get neutral losses for each monomer, add to dict if present
    loss_dict = {}
    for monomer in monomer_ids:
        potential_losses = retrieve_loss_products_monomer_id(
            monomer_id=monomer,
            polymer=polymer)
        if potential_losses:
            loss_dict[monomer] = potential_losses

    return loss_dict

def retrieve_loss_products_sequence(sequence, polymer):
    """
    Takes a sequence string, Polymer object and returns dict of monomer ids
    and associated potential neutral losses
    Args:
        sequence (str): sequence id string
        polymer (Polymer): Polymer object

    Returns:
        Dict[str, List[float]]: keys = monomer id strings, floats = list of
            neutral loss masses for monomer.
    """

    #  list of monomer ids
    monomers = sorted(return_monomer_ids_sequence(
        sequence=sequence,
        return_modified=True,
        return_set=False
    ))

    #  get neutral losses for each monomer, add to dict if present
    loss_dict = {}
    for monomer in monomers:
        potential_losses = retrieve_loss_products_monomer_id(
            monomer_id=monomer,
            polymer=polymer)
        if potential_losses:
            loss_dict[monomer] = potential_losses

    return loss_dict

def retrieve_loss_products_monomers_list(
    monomer_id_list,
    polymer
):
    """
    Retrieves neutral loss masses for monomer sidechains from list of
    monomer id strings

    Args:
        monomer_id_list (List[str]): list of monomer id strings
        polymer (Polymer): Polymer object

    Returns:
        Dict[str, List[float]]: dict of monomer id strings and associated
            lists of potential neutral loss masses
    """
    #  get neutral losses for each monomer, add to dict if present
    loss_dict = {}
    for monomer in monomer_id_list:
        potential_losses = retrieve_loss_products_monomer_id(
            monomer_id=monomer,
            polymer=polymer)
        if potential_losses:
            loss_dict[monomer] = potential_losses

    return loss_dict

@lru_cache(maxsize=20)
def retrieve_loss_products_monomer_id(monomer_id, polymer):
    """
    Takes monomer_id and Polymer object, retrieves neutral losses for Monomer
    object with matching id.

    Args:
        monomer_id (str): monomer id code. NOTE: this can be one-letter code
            for unmodified monomers or modified monomer string for sidechain
            modified monomer
        polymer (Polymer): Polymer object.

    Returns:
        List[float]: list of neutral masses corresponding to neutral losses for
            target monomer.
    """

    #  retrieve appropriate Monomer object from Polymer object
    monomer = list(filter(lambda x: x.id == monomer_id[0], polymer.monomers))[0]

    #  retrieve loss products as defined in Monomer object, if None return None
    loss_products = monomer.loss_products
    if not loss_products:
        return []

    final_neutral_losses = []

    #  iterate through losses and convert to floats for absolute mass values.
    #  Look up in FUNCTIONAL_GROUPS if losses are given as func group strings
    #  (e.g. "H2O")
    for product in loss_products:
        try:
            product = float(product)
            final_neutral_losses.append(product)
        except ValueError:
            if product[0] == "-":
                product = -FUNCTIONAL_GROUPS[product[1::]]
            else:
                product = FUNCTIONAL_GROUPS[product]
            final_neutral_losses.append(product)

    #  if monomer is unmodified, return standard neutral losses
    if len(monomer_id) == 1:
        return final_neutral_losses

    #  retrieve sidechain modification strings
    modifications = remove_substrings_sequence(
        sequence=monomer_id,
        substring_markers={'(': ')'},
        return_substrings=True
    )

    #  iterate through sidechain modifications and check if any prevent neutral
    #  losses
    for modification in modifications:
        mod_list = polymer.modifications[monomer_id[0]]
        mod_object = list(filter(
            lambda x: x.id == modification[1:-1],
            mod_list))[0]

        if mod_object.disrupt_neutral_loss:
            return []

    return final_neutral_losses

def add_sidechain_neutral_loss_products_sequence(
    sequence,
    sequence_mass,
    polymer,
    params,
    ms_level=1
):
    """
    Takes a sequence, its corresponding neutral monoisotopic mass, and returns
    a list of neutral loss products for sequence.

    Args:
        sequence (str): sequence string
        sequence_mass (float): neutral monoisotopic mass for sequence
        polymer (Polymer): Polymer object
        params (Parameters): Parameters object
        ms_level (int, optional): level of MS data sequence is to eventually be
            screened for in - this parameter is specified because different
            neutral loss fragmentation limits are placed on MS1 and MS2 data.
            Defaults to 1.

    Returns:
        List[float]: list of neutral loss products (with loss product =
            full_neutral_mass - neutral_loss).
    """
    #  init list to store products of neutral loss fragmentations
    neutral_loss_products = []

    if ms_level == 1:
        max_neutral_losses = params.silico.ms1.max_neutral_losses
    else:
        max_neutral_losses = params.silico.ms2.max_neutral_losses
    if not max_neutral_losses:
        return neutral_loss_products

    #  list monomer ids for sequence
    monomers = return_monomer_ids_sequence(
        sequence=sequence,
        return_modified=True,
        return_set=False
    )

    #  retrieve loss products from Polymer object, add each loss product
    #  to loss product list
    potential_losses = retrieve_loss_products_monomers_list(
        monomer_id_list=monomers,
        polymer=polymer
    )

    if not potential_losses:
        return []

    neutral_loss_products, total_losses = [], 0

    max_neutral_losses = min(max_neutral_losses, sum([
        monomers.count(monomer) for monomer in potential_losses
    ]))

    #  iterate through loss product-prone monomers and add corresponding loss
    #  products up until the maximum neutral loss limit has been reached
    for loss_monomer, loss_masses in potential_losses.items():
        occurence = monomers.count(loss_monomer)
        if total_losses < max_neutral_losses:
            for mass in loss_masses:
                for i in range(occurence):
                    if total_losses < max_neutral_losses:
                        neutral_loss_products.append(
                            sequence_mass - (mass * (i + 1)))
                        total_losses += 1

    #  return neutral loss products sorted by descending mass
    return sorted(neutral_loss_products, reverse=True)

def generate_all_sequences(
    polymer,
    params,
    sequencing=True
):
    """
    Takes a Polymer and Parameters object, returns a full list of potential
    sequences or compositions

    Args:
        polymer (Polymer): Polymer object
        params (Parameters): Parameters object
        sequencing (bool, optional): Specifies whether to return list of unique
            sequence or composition strings; if True, returns sequences, if
            False, returns compositions. Defaults to True.

    Returns:
        List[str]: list of sequence or composition strings.
    """

    #  if isomeric target sequences are specified in input parameters, return
    #  list of all sequences that are isomeric to one or more target sequence
    if params.silico.isomeric_targets:
        return get_isomeric_seqs(
            target_sequences=params.silico.isomeric_targets,
            sequencing=sequencing
        )

    #  return list of all possible sequences or compositions (if
    #  sequencing=False) that fit constraints of input parameters
    return get_sequences(
        polymer=polymer,
        params=params,
        sequencing=sequencing
    )

def get_sequences(
    polymer,
    params,
    sequencing
):
    """
    Get list of unique sequence strings from Polymer and Parameters objects.
    Args:
        polymer (Polymer): Polymer object
        params (Parameters): Parameters object
        sequencing (bool): Specifies whether to return list of unique
            sequence or composition strings; if True, returns sequences, if
            False, returns compositions.

    Returns:
        List[str]: list of unique sequence strings
    """
    sequences = []

    #  get monomer id strings
    monomer_ids = [
        monomer.id for monomer in polymer.monomers
    ]

    #  iterate through sequence length range, add every sequence permutation
    #  of each length to sequences list
    for length in range(
            params.silico.min_length,
            params.silico.max_length + 1):
        sequences.extend([
            ''.join(combination)
            for combination in
            list(itertools.product(monomer_ids, repeat=length))
        ])

    if not sequencing:
        return list(set(
            ["".join(sorted(seq)) for seq in sequences]))
    return list(set(sequences))

def get_compositions(
    polymer,
    params
):
    """Used by generate_all_sequences."""
    compositions = []

    #  get monomer id strings
    monomer_ids = [
        monomer.id for monomer in polymer.monomers
    ]

    for length in range(
            params.silico.min_length,
            params.silico.max_length + 1):
        compositions.extend([
            ''.join(sorted(combination))
            for combination in
            itertools.combinations_with_replacement(monomer_ids, length)
        ])

    #  sort compositions and return
    return list(set([''.join(sorted(item)) for item in compositions]))

def get_isomeric_seqs(target_sequences, sequencing=True):
    """
    Generates and returns lists of all sequences that are isomeric to one or
    more target in target_sequences.

    Args:
        target_sequences (List[str]): list of sequence strings for target
            sequences.
        sequencing (bool, optional): specifies whether to return composition
            strings or unique sequence strings. Defaults to True.

    Returns:
        List[str]: list of all sequence string permutations that are isomeric
            to one or more isomeric targets.
    """

    #  init list to store isomeric sequences
    isomeric_sequences = []

    #  iterate through isomeric targets, retrieve monomer ids and all possible
    #  sequence permutations
    for target in target_sequences:
        monomers = return_monomer_ids_sequence(
            sequence=target,
            return_modified=True,
            return_set=False
        )

        #  if returning unique sequences, get all possible sequence permutations
        #  for target monomers
        if sequencing:
            isomeric_sequences.extend([
                "".join(x) for x in itertools.permutations(monomers)
            ])

        #  returning compositions, so treat isomers as identical composition
        #  strings
        else:
            isomeric_sequences.append("".join(sorted(monomers)))

    #  return list of sequences isomeric to one or more target, sorted by
    #  sequence length
    return list(set(sorted(
        isomeric_sequences,
        key=lambda x: len(x)
    )))