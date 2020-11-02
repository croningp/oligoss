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
        monomer_id (str): monomer id string, including any sidechain or terminal
            modifications
        polymer (Polymer): instance of Polymer class

    Returns:
         float: neutral monoisotopic mass of monomer
    """

    # get monomer mass if terminal modification
    if monomer_id[0] == '[':
        monomer = list(
            filter(lambda x: x.id == monomer_id[5], polymer.monomers))[0]

    #  get appropriate Monomer object from Polymer object
    else:
        monomer = list(
            filter(lambda x: x.id == monomer_id[0], polymer.monomers))[0]

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
    associated sidechain and/or terminal modifications).

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
    terminal_modifications = {}

    if '[' in sequence:

        # create dictionary of removed modifications to add back later
        terminal_modifications = return_terminal_mods(sequence=sequence)
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

    # if terminal modifications were removed, add them back into the list
    if terminal_modifications:
        for terminus, mod in terminal_modifications.items():
            if terminus == 0:
                monomers[0] = f'[{mod}]{monomers[terminus]}'
            if terminus == -1:
                monomers[-1] = f'{monomers[terminus]}[{mod}]'

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

def return_terminal_mods(sequence):
    """ This function creates a dictionary listing the modifications present
    at terminal monomers for a given sequence. Used in get_modified_monomer_ids.

    Args:
        sequence (str): terminally modified sequence string.

    Returns:
        terminal_modifications (Dict[int, str]): key is terminus position and
            value is three letter modification string.
    """

    terminal_modifications = {}

    if sequence[0] == '[':
        terminal_modifications[0] = str(sequence[1:4])

    if sequence[-1] == ']':
        terminal_modifications[-1] = str(sequence[-4:-1])

    return terminal_modifications

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
    #  get adducts for set ms_level
    if ms_level == 1:
        min_z, max_z = params.silico.ms1.min_z, params.silico.ms1.max_z
        adducts = params.silico.ms1.adducts
    else:
        min_z, max_z = params.silico.ms2.min_z, params.silico.ms2.max_z
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
            for unmodified monomers or modified monomer string for sidechain or
            terminally modified monomers
        polymer (Polymer): Polymer object.

    Returns:
        List[float]: list of neutral masses corresponding to neutral losses for
            target monomer.
    """
    # get monomer mass if terminal modification
    if monomer_id[0] == '[':
        monomer = list(
            filter(lambda x: x.id == monomer_id[5], polymer.monomers))[0]

    else:
        #  retrieve appropriate Monomer object from Polymer object
        monomer = list(
            filter(lambda x: x.id == monomer_id[0], polymer.monomers))[0]

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

        # account for sequences modified at 0 terminus 0
        if monomer_id[0] == '[':
            mod_mon = monomer_id[5]
        else:
            mod_mon = monomer_id[0]

        mod_list = polymer.modifications[mod_mon]
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

    if ms_level == 1:
        max_neutral_losses = params.silico.ms1.max_neutral_losses
    else:
        max_neutral_losses = params.silico.ms2.max_neutral_losses
    if max_neutral_losses == 0:
        return []

    #  list monomer ids for sequence
    monomers = return_monomer_ids_sequence(
        sequence=sequence,
        return_modified=True,
        return_set=False
    )

    #  retrieve loss products from Polymer object, add each loss product
    #  to loss product list
    potential_unique_losses = retrieve_loss_products_monomers_list(
        monomer_id_list=monomers,
        polymer=polymer
    )

    #  no loss product-prone sidechains so sashay away
    if not potential_unique_losses:
        return []

    #  if no cap on max_neutral_losses has been set, limit this only by number
    #  possible due to intrinsic properties of sequence monomers
    if not max_neutral_losses:
        max_neutral_losses = len(sequence)**2

    #  init list to store loss products, counter for losses
    neutral_loss_products, total_losses = [sequence_mass], 0

    #  iterate through loss product-prone monomers and add corresponding loss
    #  products up until the maximum neutral loss limit has been reached
    for loss_monomer, loss_masses in potential_unique_losses.items():
        occurence = monomers.count(loss_monomer)
        monomer_losses = []
        for loss_mass in loss_masses:
            for i in range(occurence):
                final_loss = loss_mass * (i + 1)
                monomer_losses.extend([
                    mass - final_loss
                    for mass in neutral_loss_products
                ])
        for loss_product in monomer_losses:
            if total_losses < max_neutral_losses:
                if loss_product not in neutral_loss_products:
                    neutral_loss_products.append(loss_product)
                    total_losses += 1

    #  remove intact sequence mass from loss products list
    neutral_loss_products.remove(sequence_mass)

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
            polymer=polymer,
            sequencing=sequencing
        )

    #  generate list of all possible sequences or compositions (if
    #  sequencing=False) that fit constraints of input parameters
    compositions = get_sequences(
        polymer=polymer,
        params=params,
        sequencing=sequencing)

    silico = params.silico.ms1

    if params.silico.modifications:
        compositions = add_modifications_strings(
            sequence_list=compositions,
            modifications=params.silico.modifications,
            universal_terminal_mod=silico.universal_terminal_modifications)
    return (comp for comp in compositions)

def get_sequences(
    polymer,
    params,
    sequencing
):
    """
    Get unique sequence strings from Polymer and Parameters objects.
    Args:
        polymer (Polymer): Polymer object
        params (Parameters): Parameters object
        sequencing (bool): Specifies whether to return list of unique
            sequence or composition strings; if True, returns sequences, if
            False, returns compositions.

    Yields:
        str: sequence string (if sequencing == True) or composition string (if
            sequencing == False).
    """
    #  init list to store sequence generators. NOTE: this will be populated with
    #  list of Generator[str] type objects, not individual sequence strings.
    sequences = []

    #  get monomer id strings
    monomer_ids = [
        monomer.id for monomer in polymer.monomers
    ]

    #  iterate through sequence length range, append Generator[str] objects for
    #  sequences at each length to sequences list.
    for length in range(
            params.silico.min_length,
            params.silico.max_length + 1):
        sequences.append(get_seqs_fixed_len(
            length=length, monomer_ids=monomer_ids))

    #  in case of duplicates, keep track of which sequence strings have been
    #  yielded. This has an obvious memory cost but will avoid the disaster of
    #  duplicates. NOTE: this is only required when sequencing=False as it is
    #  only compositions that can potentially be duplicated.
    final_seqs = []
    for seq_generator in sequences:
        for seq_str in seq_generator:
            if sequencing:
                yield seq_str
            else:
                composition = "".join(sorted(seq_str))
                if composition not in final_seqs:
                    final_seqs.append(composition)
                    yield composition

def get_seqs_fixed_len(length, monomer_ids, isomeric_target=None):
    """
    Takes monomer ids and target length, returns generator of all possible
    sequence permutations of target length with monomer_ids.

    Args:
        length (int): target length of sequence
        monomer_ids (List[str]): list of monomer id strings.
        isomeric_targets (Optional[str]): target composition. If specified, only
            sequences isomeric to target will be yielded. Defaults to None.
    Yields:
        str: sequence string.
    """
    for seq in itertools.product(monomer_ids, repeat=length):
        if not isomeric_target:
            yield "".join(seq)
        else:
            if "".join(sorted(seq)) == isomeric_target:
                yield "".join(seq)

def get_isomeric_seqs(target_sequences, polymer, sequencing=True):
    """
    Generates and returns lists of all sequences that are isomeric to one or
    more target in target_sequences.

    Args:
        target_sequences (List[str]): list of sequence or composition strings
            for targetting isomers.
        polymer (Polymer): Polymer object.
        sequencing (bool, Optional): specifies whether to yield unique sequence
            strings (True) or just compositions (False). Defaults to False.

    Yields:
        str: sequence or composition string permutations that are
            isomeric to one or more isomeric targets.
    """

    #  init list to store isomeric sequences
    #  NOTE: this will be a list of generators, not individual sequences
    isomeric_sequences = []

    #  iterate through isomeric targets, retrieve monomer ids and all possible
    #  sequence permutations
    target_sequences = set(
        ["".join(sorted(target)) for target in target_sequences]
    )

    for target in target_sequences:
        monomer_ids = return_monomer_ids_sequence(
            sequence=target,
            return_modified=True,
            return_set=False
        )

        isomeric_sequences.append(get_seqs_fixed_len(
            length=len(monomer_ids),
            monomer_ids=list(set(monomer_ids)),
            isomeric_target=target
        ))

    #  init list to store seqs that have been yielded. NOTE: this is only used
    #  if sequencing == False as it is to remove the possibility of repeat
    #  compositions being returned. There will be a memory hit, but not for
    #  the large sequence pools where sequencing == True
    previous_seqs = []

    #  iterate through sequence generators
    for seq_generator in itertools.chain(tuple(isomeric_sequences)):
        #  iterate through sequences in generator and yield sequence strings
        for seq in seq_generator:
            if not sequencing:
                comp = "".join(sorted(seq))
                if comp not in previous_seqs:
                    previous_seqs.append(comp)
                    yield comp
            else:
                yield seq

def generate_isomers_target(target, polymer):
    """
    Takes a target sequence and returns generator of all sequences isomeric to
    target.

    Args:
        target (str): target sequence string.
        polymer (Polymer): Polymer object.

    Yields:
        str: sequence string
    """

    #  get list of monomer ids, including repeated monomers
    # e.g. for target "AK(Ole)AAV" monomers = ["A", "K(Ole)", "A", "A", "V"]
    monomers = return_monomer_ids_sequence(
        sequence=target,
        return_modified=True,
        return_set=False)

    #  if returning unique sequences, get all possible sequence permutations
    #  for target monomers
    for sequence in itertools.permutations(monomers):
        yield "".join(sequence)

def add_modifications_strings(
    sequence_list,
    modifications,
    universal_terminal_mod=True
):
    """ This function adds terminal modifications to sequences in the
    sequence list. If applicable, the sequence list will already include
    sidechain modifications at this point.

    Args:
        sequence_list (iterator): iterator (list or generator) of sequences.
        modifications (dict): dictionary of all modification targets and their
            modifications (sidechain and terminal).
        universal_terminal_mod (bool, optional): Whether or not the terminal
            modification is always present in sequences. Defaults to True.

    Returns:
        Generator[str]: generator of modified sequences.
    """

    modified_sequences = sequence_list

    terminal_mods = {
        k: v for k, v in modifications.items() if k in ['-1', '0']}

    if not terminal_mods:
        return sequence_list

    if terminal_mods:
        modified_sequences = add_terminal_mods_strings(
            sequence_list=modified_sequences,
            terminal_modifications=terminal_mods,
            universal_mod=universal_terminal_mod)

    return list(set(modified_sequences))

def add_terminal_mods_strings(
    sequence_list, terminal_modifications, universal_mod
):
    """ This function takes a list of sequence strings and adds terminal
    modifications of the format [mod] to the desired terminus.

    Args:
        sequence_list (List[str]): list of sequence strings to be modified.
        terminal_modifications (Dict[str, List[str]]): keys are target termini
            and their corresponding key is a list of possible terminal
            modifications.
        universal_mod (Bool): states whether or not the terminal modification is
            universal / all sequences must be terminally modified.

    Returns:
        List[str]: list of terminally modified sequences.
    """

    if not terminal_modifications:
        return sequence_list

    terminal_modified = sequence_list

    for target, mods in terminal_modifications.items():

        if type(mods) != list:
            mods = [mods]
        single_target = []

        for sequence in terminal_modified:

            for mod in mods:

                # add terminally modified string to list if not modified already
                if target == '0' and sequence[0] != '[':
                    single_target.append(f'[{str(mod)}]{sequence}')

                if target == '-1' and sequence[-1] != ']':
                    single_target.append(f'{sequence}[{str(mod)}]')

        terminal_modified.extend(single_target)

    if universal_mod:
        return list(
            set([t for t in terminal_modified if t.count('[') == len(target)]))

    return list(set(terminal_modified))

def get_composition(sequence):
    """ This function generates an alphabetically ordered composition string
    from a sequence string, including those with terminal and/or sidechain
    modifications.

    Args:
        sequence (str): sequence string.

    Returns:
        str: compositional sequence string.
    """
    if '[' in sequence:
        modified_mons = {}
        terminal_mods = return_terminal_mods(sequence=sequence)
        seq_list = get_modified_monomer_ids(sequence=sequence, return_set=False)

        for terminus in terminal_mods:

            # remove terminal modifications but store their original position
            modified_mons[terminus] = seq_list[terminus]
            del seq_list[terminus]

        # sort un-terminally-modified sequence
        sorted_seq = sorted(seq_list)

        # add terminal modification back into sequence
        for terminus, mod in modified_mons.items():
            if terminus == 0:
                sorted_seq.insert(0, mod)
            if terminus == -1:
                sorted_seq.append(mod)
        return ''.join(sorted_seq)

    return ''.join(sorted(
        get_modified_monomer_ids(sequence=sequence, return_set=False)))
