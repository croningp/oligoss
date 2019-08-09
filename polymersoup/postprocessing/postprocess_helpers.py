"""
This file contains core operations required by postprocessing module
"""

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
    core fragment series

    Args:
        confirmed_fragments (list): list of confirmed fragment ids
        core_fragment_series (list): list of fragment one letter codes for
                        series that are to be used in calculation
        optional_fragments (list): list of specific fragment ids to
                        exclude from calculation ONLY IF they are not in the
                        list of confirmed fragments. Include in the calculation
                        if they are in the list
        excluded_fragments (list): list of specific fragment ids to
                        exclude from subsequence coverage in ALL cases
    Returns:
        final_sequence_coverage (float): average % sequence coverage for all #
                        core fragment series

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
        silico_core_fragments = [
            frag for frag in silico_fragments
            if frag[0] == core_series
        ]

        # get confirmed fragments that belong to current core series
        core_confirmed = [
            frag for frag in confirmed_fragments
            if frag[0] == core_series
        ]

        # get optional fragments that belong to current core series
        optional_core_fragments = [
            frag for frag in silico_core_fragments
            if frag in optional_fragments
        ]

        # get list of core fragment positions from in silico fragments,
        # sorted in ascending order (e.g. ['y1', 'y3'] becomes [1, 3])
        core_indeces = sorted([
            int(frag[1::])
            for frag in silico_core_fragments
            if frag not in excluded_fragments
        ])

        # get list of confirmed fragment positions from confirmed fragments,
        # sorted in ascending order (e.g. ['y1', 'y3'] becomes [1, 3])
        confirmed_indeces = sorted([
            int(frag[1::])
            for frag in confirmed_fragments
            if frag not in excluded_fragments
        ])

        # calculate largest POSSIBLE continuous block of fragments for this
        # core series
        max_sequence_coverage = calculate_subsequence_coverage(
            core_indeces
        )

        # calculate largest OBSERVED continuous block of fragmentsfor this
        # core series
        obs_sequence_coverage = calculate_subsequence_coverage(
            confirmed_indeces
        )

        # calculate percentage of sequence covered by confirmed fragments in
        # this core series
        sequence_coverage = (obs_sequence_coverage/max_sequence_coverage)*100

        # add sequence coverage of this core series to coverage_dict
        coverage_dict[core_series] = sequence_coverage

    # final sequence coverage is an average of the % coverage of all core
    # fragment series
    final_sequence_coverage = sum(coverage_dict.values())/len(coverage_dict)

    # return final sequence coverage (units = %)
    return final_sequence_coverage

def calculate_subsequence_coverage(
    fragment_indeces
):
    """
    Takes a list of fragment indeces and returns sequence coverage (i.e.
    number of fragments in longest continuous block of fragments)

    Args:
        fragment_indeces (list of ints): list of fragment positions in a given
            series (ints)

    Returns:
        sequence_coverage (int): number of fragments in longest continuous
            block of fragments
    """
    # set sequence coverage to 0 before iterating through fragment indeces
    sequence_coverage = 0

    # iterate through fragment indeces, performing a walk through the list of
    # indeces until there is a break in the series of indeces
    for i in range(0, len(fragment_indeces)):

        # choose starting index
        start_fragment = fragment_indeces[i]

        # initialse coverage from start_fragment to
        sub_sequence_length = 0

        # increment sub_sequence_length as long as there the fragment series
        # is uninterrupted
        while start_fragment + 1 in fragment_indeces:
            sub_sequence_length += 1
            start_fragment += 1

        # reset sequence coverage if longest block of continuous fragment
        # indeces exceeds previous value
        if sub_sequence_length > sequence_coverage:
            sequence_coverage = sub_sequence_length

    # return sequence_coverage (i.e. number of fragments in longest continuous
    # block)
    return sequence_coverage

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
    """
    This function takes in silico and confirmed fragments for a target
    sequence and assigns a confidence score using the confidence scoring
    criteria set

    Args:
        insilico_fragments (list of strings): list of fragment ids that are
                    theoretically possible for sequence
        confirmed_fragments (list of strings): list of fragment ids for
                    fragments that have been confirmed from mass spec data
        core_fragment_series (list of strings): list of fragment series one
                    letter codes for fragment series used in calculating %
                    confirmed fragments and / or fragment sequence coverage
        optional_fragments (list of strings): list of fragment ids for
                    fragments that can be excluded from consideration IF they
                    are not confirmed
        exclude_fragments (list of strings): list of fragment ids for
                    fragments that are to be excludeed from consideration
                    whether confirmed or not. These will never be used in
                    confidence calculation
        essential_fragments (list of strings): list of fragment ids for
                    fragments whose presence is required for confidence
                    assignment to exceed a threshold (essential_fragment_cap)
        essential_fragment_cap (float): upper limit on assigned confidence
                    threshold for a sequence if any essential fragments are
                    missing
        sequence_coverage_weight (float): weighting assigned to sequence
                    coverage for confidence calculation (given as a decimal
                    fraction - i.e. MUST BE BETWEEN 0 AND 1)

    Returns:
        confidence (float): final confidence assignment (in % )
    """

    # if any essential fragments are missing, restrict confidence assignment
    # to an upper limit of essential_fragment_cap
    missing_essentials = []

    if essential_fragments:
        missing_essentials.extend([
            frag for frag in essential_fragments
            if frag not in confirmed_fragments
        ])

    # default cap on confidence assignment is 100%; only reduce this if
    # essential fragments are missing
    confidence_cap = 100
    if missing_essentials:
        confidence_cap = essential_fragment_cap

    # add any optional fragments that have not been confirmed to list of
    # fragments to exclude from consideration
    exclude_fragments.extend([
        frag for frag in optional_fragments
        if frag not in confirmed_fragments
    ])

    # remove exclude_fragments from in silico fragments
    insilico_fragments = [
        frag for frag in insilico_fragments
        if frag not in exclude_fragments
    ]

    # remove exclude_fragments from confirmed fragments
    confirmed_fragments = [
        frag for frag in confirmed_fragments
        if frag not in exclude_fragments
    ]

    # calculate percentage of theoretical (in silico) fragments that have been
    # observed in confirmed_fragments
    n_confirmed, n_silico = len(confirmed_fragments), len(insilico_fragments)
    percentage_found_fragments = (n_confirmed / n_silico)*100

    # set sequence coverage to 0; this will only be calculated if
    # sequence_coverage_weight is greater than 0
    sequence_coverage = 0

    # check whether sequence coverage is to be taken into account; if so,
    # calculate sequence coverage
    if sequence_coverage_weight > 0:
        sequence_coverage = get_subsequence_coverage(
            insilico_fragments,
            confirmed_fragments,
            core_fragment_series,
            optional_fragments,
            exclude_fragments
        )

    # calculate final confidence by combining WEIGHTED values for percentage
    # confirmed fragments and sequence coverage
    confidence = percentage_found_fragments*(1-sequence_coverage_weight)
    confidence += sequence_coverage*sequence_coverage_weight

    # return final confidence score for sequence, making sure to apply
    # confidence cap if any essential fragments are missing
    return min([confidence, confidence_cap])

def trim_EIC(
    EIC,
    min_Rt,
    max_Rt,
    min_absolute_intensity,
    min_relative_intensity
):
    """
    This function trims an EIC, removing data points that violate the
    thresholds specified

    Args:
        EIC (list): list of retention time + intensity sub_lists in format:
            [[Rt, I]...] where Rt = retention time, I = signal intensity at
            Rt
        min_Rt (float): minimum retention time, anything below which is removed
        max_Rt (float): maximum retention time, anything above which is removed
        min_absolute_intensity (float): minimum intensity of data points,
            anything below which is removed
        min_relative_intensity (float): minimum RELATIVE intensity of data
            points as a decimal fraction of most intense point in the EIC


    Returns:
        EIC: trimmed EIC in same format as input, but with data trimmed
    """

    # get intensity of most intense peak in the EIC
    max_intensity = sorted(lambda Rt_I: Rt_I[1])[1]

    # if max_Rt is not defined, set max_Rt to latest retention time point in
    # the EIC
    if not max_Rt:
        max_Rt = sorted(lambda Rt_I: Rt_I[0], EIC)[-1][0]

    # if min_Rt is not defined, set min_Rt to earliest retention time point in
    # the EIC
    if not min_Rt:
        min_Rt = sorted(lambda Rt_I: Rt_I[0], EIC)[0][0]

    # remove data that fall outwith retention time range
    EIC = filter(
        lambda Rt_I: Rt_I[0] >= min_Rt and Rt_I[0] <= max_Rt, EIC
    )

    # remove data below minimum absolute intensity, if this is specified
    if min_absolute_intensity:
        EIC = filter(
            lambda Rt_I: Rt_I[1] >= min_absolute_intensity, EIC
        )

    # remove data below minimum relative intensity (i.e. % intensity of most
    # intense peak - expressed as a decimal fraction in min_relative_intensity
    # variable)
    if min_relative_intensity:
        EIC = filter(
            lambda Rt_I: Rt_I[1] >= min_relative_intensity*max_intensity, EIC
        )

    # return trimmed EIC
    return EIC

def get_Rt_I_from_EIC(
    EIC,
    n_targets=1,
    min_Rt=None,
    max_Rt=None,
    min_absolute_intensity=None,
    min_relative_intensity=None,
    Rt_bin=0.2
):

    EIC = trim_EIC(
        EIC=EIC,
        min_Rt=min_Rt,
        max_Rt=max_Rt,
        min_absolute_intensity=min_absolute_intensity,
        min_relative_intensity=min_relative_intensity
    )
