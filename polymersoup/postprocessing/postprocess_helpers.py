def assign_confidence_sequences(
    silico_dict,
    confirmed_dict,
    postprocess_params,
    ssw
):
    """ This function assigns confidence scores to each sequence in the
    confirmed dict and returns a dictionary of sequences and their confidence
    scores.

    Arguments:
        silico_dict (dict) -- dictionary of all possible sequences
        (silico dictionary).
        confirmed_dict (dict) -- dictionary of confirmed sequences.
        postprocess_params (dict) -- dictionary of postprocessing parameters,
        retrieved from input parameters JSON file.

    Returns:
        confidence_dict (dict) -- dictionary of sequences and their confidence
        scores.
    """
    # remove sequences from in silico dict that do not have any confirmed
    # fragments
    silico_dict = {
        sequence: subdict for sequence, subdict in silico_dict.items()
        if sequence in confirmed_dict.keys()}

    # initiate dict to store confidence scores for sequences
    confidence_dict = {}

    # iterate through sequences that have confirmed fragments, assigning each
    # a confidence score
    for sequence, confirmed_fragments in confirmed_dict.items():

        print(f"assigning confidence for {sequence}")

        # assign sequence a confidence score, using constraints specified in
        # input postprocess parameters from input parameters .json file
        confidence_assignment = assign_confidence_sequence(
            insilico_fragments=silico_dict[sequence]["MS2"],
            confirmed_fragments=confirmed_fragments["MS2"],
            core_fragment_series=postprocess_params.core_linear_series,
            optional_fragments=postprocess_params.optional_core_fragments,
            exclude_fragments=postprocess_params.exclude_fragments,
            essential_fragments=postprocess_params.essential_fragments,
            essential_fragment_cap=postprocess_params.dominant_signature_cap,
            sequence_coverage_weight=ssw)

        # add sequence and its final confidence score to confidence_dict
        confidence_dict[sequence] = confidence_assignment

    # return dictionary of sequences and their confidence scores
    return confidence_dict

def assign_confidence_sequence(
    insilico_fragments,
    confirmed_fragments,
    core_fragment_series,
    optional_fragments,
    exclude_fragments,
    essential_fragments,
    essential_fragment_cap,
    sequence_coverage_weight
):
    """ This function takes in silico and confirmed fragments for a target
    sequence and assigns a confidence score using the confidence scoring
    criteria set.

    Args:
        insilico_fragments (dict[str]): dict of fragment ids and their
        corresponding masses that are theoretically possible for sequence.
        confirmed_fragments (list[str]): list of fragment ids for fragments that
        have been confirmed from mass spec data.
        core_fragment_series (list[str]): list of fragment series one letter
        codes for fragment series used in calculating % confirmed fragments and
        / or fragment sequence coverage.
        optional_fragments (list[str]): list of fragment ids for fragments that
        can be excluded from consideration IF they are not confirmed.
        exclude_fragments (list[str]): list of fragment ids for fragments that
        are to be excludeed from consideration whether confirmed or not.
        These will never be used in confidence calculation.
        essential_fragments (list[str]): list of fragment ids for fragments
        whose presence is required for confidence assignment to exceed a
        threshold (essential_fragment_cap).
        essential_fragment_cap (float): upper limit on assigned confidence
        threshold for a sequence if any essential fragments are missing.
        sequence_coverage_weight (float): weighting assigned to sequence
        coverage for confidence calculation (given as a decimal fraction - i.e.
        MUST BE BETWEEN 0 AND 1).

    Returns:
        confidence (float): final confidence assignment (in %).
    """
    insilico_fragment_dict = {
        'core': [k for k in insilico_fragments.keys() if k != 'signatures'],
        'signatures': [
            s for s in insilico_fragments['signatures']
            if s != 'terminal_modifications'],
        'terminal_modifications': [t for t in insilico_fragments[
            'signatures']['terminal_modifications'].keys()]}

    confirmed_fragment_dict = {
        'core': [k for k in confirmed_fragments['core'].keys()],
        'signatures': [
            s for s in confirmed_fragments['signatures']
            if s != 'terminal_modifications'],
        'terminal_modifications': [t for t in confirmed_fragments[
            'signatures']['terminal_modifications'].keys()]}

    # work out percentage found fragments from core series
    percentage_found_fragments = get_percent_found_fragments(
        insilico_fragments=insilico_fragment_dict,
        confirmed_fragments=confirmed_fragment_dict,
        core_fragment_series=core_fragment_series,
        optional_fragments=optional_fragments,
        exclude_fragments=exclude_fragments)

    # initiate list to store missing essential fragments
    missing_essentials = []

    # check that all essential fragments are in insilico fragments
    if essential_fragments:
        essential_fragments = list(
            filter(
                lambda frag: frag in insilico_fragments, essential_fragments))

    # if any essential fragments are missing, restrict confidence assignment
    # to an upper limit of essential_fragment_cap
        missing_essentials = [frag for frag in essential_fragments
                              if frag not in confirmed_fragments]

    # default cap on confidence assignment is 100%; only reduce this if
    # essential fragments are missing
    confidence_cap = 100
    if missing_essentials:
        confidence_cap = essential_fragment_cap

    # set sequence coverage to 0; this will only be calculated if
    # sequence_coverage_weight is greater than 0
    sequence_coverage = 0

    # check whether sequence coverage is to be taken into account; if so,
    # calculate sequence coverage
    if sequence_coverage_weight > 0:
        sequence_coverage = get_subsequence_coverage(
            silico_fragments=insilico_fragment_dict,
            confirmed_fragments=confirmed_fragment_dict,
            core_fragment_series=core_fragment_series,
            optional_fragments=optional_fragments,
            excluded_fragments=exclude_fragments)

    # calculate final confidence by combining WEIGHTED values for percentage
    # confirmed fragments and sequence coverage
    confidence = percentage_found_fragments * (1 - sequence_coverage_weight)
    confidence += sequence_coverage * sequence_coverage_weight

    # return final confidence score for sequence, making sure to apply
    # confidence cap if any essential fragments are missing
    return min([confidence, confidence_cap])

def get_percent_found_fragments(
    insilico_fragments,
    confirmed_fragments,
    core_fragment_series,
    optional_fragments,
    exclude_fragments
):
    """ Takes a list of in silico fragments for a sequence, confirmed fragments
    and one letter codes for core fragments, and returns % core fragments that
    have been confirmed for sequence.

    Args:
        insilico_fragments (dict{list}): dict of fragment types and a list of
        fragment ids for in silico (i.e.theoretical) sequence fragment.
        confirmed_fragments (list): list of fragment ids for confirmed (i.e.
        observed) sequence fragments.
        core_fragment_series (list): list of ONE LETTER CODES for fragment
        series that are to be used to calculate percentage found fragments.
        optional_fragments (list): list of specific fragment ids to exclude
        from calculation IF they have not been confirmed.
        exclude_fragments (list): list of specific fragment ids to exclude
        from calculation UNDER ANY CIRCUMSTANCES, including if they have
        been confirmed.

    Returns:
        float: % fragments that have been confirmed for a sequence.
    """
    # get list of in silico fragments that belong to core fragment series
    core_silico = insilico_fragments["core"]

    # get list of confirmed fragments that belong to core fragment series
    core_confirmed = confirmed_fragments["core"]

    # check for fragments to exclude from confidence calculation; if any are
    # specified, remove these fragments from silico and confirmed fragment
    # lists if they are present
    if exclude_fragments:
        core_silico = list(filter(
            lambda frag: frag not in exclude_fragments, core_silico))

        core_confirmed = list(filter(
            lambda frag: frag not in exclude_fragments, core_confirmed))

    # check for optional fragments; if specified, remove any optional fragments
    # that have not been confirmed from silico fragment list
    if optional_fragments:
        missing_optional_fragments = [frag for frag in core_silico if (
            frag not in core_confirmed and frag in optional_fragments)]

        core_silico = list(filter(
            lambda frag: frag not in missing_optional_fragments, core_silico))

    print(f"core_silico = {core_silico}")
    print(f"core_confirmed = {core_confirmed}")

    # count number of theoretical (silico) and confirmed fragments
    n_silico, n_confirmed = len(core_silico), len(core_confirmed)

    print(f"n_confirmed = {n_confirmed}")
    print(f"n_silico = {n_silico}")

    # return 0 if no fragments are confirmed
    if n_confirmed == 0:
        return 0

    # return % found fragments to 2 decimal places
    return round((n_confirmed / n_silico) * 100, 2)

def get_subsequence_coverage(
    silico_fragments,
    confirmed_fragments,
    core_fragment_series,
    optional_fragments,
    excluded_fragments
):
    """
    This function calculates the % a sequence that is 'covered' by continuous
    blocks of core fragment series, returning an average % coverage of all
    core fragment series.

    Args:
        confirmed_fragments (list): list of confirmed fragment ids.
        core_fragment_series (list): list of fragment one letter codes for
        series that are to be used in calculation.
        optional_fragments (list): list of specific fragment ids to exclude from
        calculation ONLY IF they are not in the list of confirmed fragments.
        Include in the calculation if they are in the list.
        excluded_fragments (list): list of specific fragment ids to exclude from
        subsequence coverage in ALL cases.
    Returns:
        final_sequence_coverage (float): average % sequence coverage for all
        core fragment series.

    """
    # initiate dict to store sequence coverages for each core fragment series;
    # dict will be in format:
    #                           {f'{core_series}' : % coverage}
    # where core_series = one letter code of fragment series, % coverage =
    # % of sequence contained within the largest continuous block of confirmed
    # fragments belonging to that series
    coverage_dict = {}

    # iterate through core fragment series, defined by their one letter codes
    for core_series in core_fragment_series:

        # get in silico fragments that belong to current core series
        silico_core_fragments = get_core_fragments(
            target_fragments=silico_fragments["core"],
            core_series=core_series)

        # get confirmed fragments that belong to current core series
        core_confirmed = get_core_fragments(
            target_fragments=confirmed_fragments["core"],
            core_series=core_series)

        # if excluded fragments, remove from consideration
        if excluded_fragments:

            core_confirmed = sorted(list(filter(
                lambda x: x not in excluded_fragments, core_confirmed)))

            silico_core_fragments = sorted(list(filter(
                lambda x: x not in excluded_fragments, silico_core_fragments)))

        # get list of core fragment positions from in silico fragments,
        # sorted in ascending order (e.g. ['y1', 'y3'] becomes [1, 3])
        core_indices = sorted(
            [int(frag[1::]) for frag in silico_core_fragments])

        # get list of confirmed fragment positions from confirmed fragments,
        # sorted in ascending order (e.g. ['y1', 'y3'] becomes [1, 3])
        confirmed_indices = sorted([int(frag[1::]) for frag in core_confirmed])

        # calculate largest POSSIBLE continuous block of fragments for this
        # core series
        max_sequence_coverage = calculate_subsequence_coverage(core_indices)

        # calculate largest OBSERVED continuous block of fragmentsfor this
        # core series
        obs_sequence_coverage = calculate_subsequence_coverage(
            confirmed_indices)

        if obs_sequence_coverage > 0:

            # calculate percentage of sequence covered by confirmed fragments in
            # this core series
            sequence_coverage = (
                obs_sequence_coverage / max_sequence_coverage) * 100
        else:
            sequence_coverage = 0

        # add sequence coverage of this core series to coverage_dict
        coverage_dict[core_series] = sequence_coverage

    # final sequence coverage is an average of the % coverage of all core
    # fragment series
    final_sequence_coverage = sum(coverage_dict.values()) / len(coverage_dict)

    # return final sequence coverage (units = %)
    return final_sequence_coverage

def get_core_fragments(
    target_fragments,
    core_series
):
    """ Takes a list of fragment ids, target core fragment series (one letter
    codes) and returns list of target fragments that fall within one or more
    of the specified core series.

    Args:
        target_fragments (list): list of fragment id strings.
        core_series (list or str): list of core fragment series one letter
        codes OR single one letter code for single core series.

    Returns:
        list: list of fragments in input target fragments that belong to one
        or more core series.
    """
    # initialise list to store core fragments
    core_fragments = []

    # make core_series a list if not already
    if type(core_series) != list:
        core_series = [core_series]

    # iterate through core series, adding fragments from target that belong
    # to one or more core series to core_fragments
    for frag_series in core_series:
        core_fragments.extend([
            frag for frag in target_fragments if frag[0] == frag_series])

    return list(set(core_fragments))

def calculate_subsequence_coverage(fragment_indices):
    """ Takes a list of fragment indeces and returns sequence coverage (i.e.
    number of fragments in longest continuous block of fragments).

    Args:
        fragment_indices (list of ints): list of fragment positions in a given
        series (ints).

    Returns:
        sequence_coverage (int): number of fragments in longest continuous
        block of fragments.
    """
    # set sequence coverage to 0 before iterating through fragment indeces
    sequence_coverage = 0

    # iterate through fragment indeces, performing a walk through the list of
    # indeces until there is a break in the series of indeces
    for i in range(0, len(fragment_indices)):

        # choose starting index
        start_fragment = fragment_indices[i]

        # initialse coverage from start_fragment to
        sub_sequence_length = 1

        # increment sub_sequence_length as long as there the fragment series
        # is uninterrupted
        while start_fragment + 1 in fragment_indices:
            sub_sequence_length += 1
            start_fragment += 1

        # reset sequence coverage if longest block of continuous fragment
        # indeces exceeds previous value
        if sub_sequence_length > sequence_coverage:
            sequence_coverage = sub_sequence_length

    # return sequence_coverage
    return sequence_coverage

def assign_isomeric_groups(sequences):
    """ This function groups isomeric sequences and assigns each group a number.

    Args:
        sequences (list[str]): list of sequences to be grouped by composition

    Returns:
        assignments (dict): Dictionary containing sequences as keys and their
        assigned isomeric group number as a value.
    """
    groups = {}

    # group sequences by composition
    for sequence in sequences:

        # temporary code - to be replaced by silico 'return monomers' function
        composition = ''.join(sorted(sequence))

        if composition in groups.keys():
            groups[composition].append(sequence)

        if composition not in groups.keys():
            groups[composition] = [sequence]

    g = 1
    assignments = {}

    # assign each sequence a isomeric group number based on composition
    for sequences in groups.values():
        for sequence in sequences:
            assignments[sequence] = g
        g += 1

    return assignments
