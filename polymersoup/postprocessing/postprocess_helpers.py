"""
This file contains core operations required by postprocessing module
"""

import pandas as pd
from pandas.io.json import json_normalize
import json
import matplotlib.pyplot as plt
import seaborn as sns


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
        silico_core_fragments = get_core_fragments(
            target_fragments=silico_fragments,
            core_series=core_series
        )

        # get confirmed fragments that belong to current core series
        core_confirmed = get_core_fragments(
            target_fragments=confirmed_fragments,
            core_series=core_series
        )
        
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
            for frag in core_confirmed
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
        sub_sequence_length = 1

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

def get_core_fragments(
    target_fragments,
    core_series
):
    """
    Takes a list of fragment ids, target core fragment series (one letter 
    codes) and returns list of target fragments that fall within one or more
    of the specified core series
    
    Args:
        target_fragments (list): list of fragment id strings
        core_series (list or str): list of core fragment series one letter 
            codes OR single one letter code for single core series 
    
    Returns:
        list: list of fragments in input target fragments that belong to one
            or more core series 
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
            frag for frag in target_fragments 
            if frag[0] == frag_series
        ])

    return list(set(core_fragments))

def get_percent_found_fragments(
    insilico_fragments,
    confirmed_fragments,
    core_fragment_series,
    optional_fragments,
    exclude_fragments
):
    """
    Takes a list of in silico fragments for a sequence, confirmed fragments 
    and one letter codes for core fragments, and returns % core fragments that
    have been confirmed for sequence 
    
    Args:
        insilico_fragments (list): list of fragment ids for in silico (i.e.
            theoretical) sequence fragment
        confirmed_fragments (list): list of fragment ids for confirmed (i.e.
            observed) sequence fragments
        core_fragment_series (list): list of ONE LETTER CODES for fragment 
            series that are to be used to calculate percentage found fragments
        optional_fragments (list): list of specific fragment ids to exclude
            from calculation IF they have not been confirmed
        exclude_fragments (list): list of specific fragment ids to exclude 
            from calculation UNDER ANY CIRCUMSTANCES, including if they have
            been confirmed 
    
    Returns:
        float: % fragments that have been confirmed for a sequence 
    """
    # get list of in silico fragments that belong to core fragment series 
    core_silico = get_core_fragments(
        target_fragments=insilico_fragments,
        core_series=core_fragment_series
    )

    # get list of confirmed fragments that belong to core fragment series 
    core_confirmed = get_core_fragments(
        target_fragments=confirmed_fragments,
        core_series=core_fragment_series
    )
    
    # check for fragments to exclude from confidence calculation; if any are
    # specified, remove these fragments from silico and confirmed fragment
    # lists if they are present 
    if exclude_fragments:
        core_silico = list(
            filter(
                lambda frag: frag not in exclude_fragments, 
                core_silico
            )
        )

        core_confirmed = list(
            filter(
                lambda frag: frag not in exclude_fragments,
                core_confirmed
            )
        )

    # check for optional fragments; if specified, remove any optional fragments
    # that have not been confirmed from silico fragment list 
    if optional_fragments:
        missing_optional_fragments = [
            frag for frag in core_silico
            if (frag not in core_confirmed 
            and frag in optional_fragments)
        ]

        core_silico = list(
            filter(
                lambda frag: frag not in missing_optional_fragments,
                core_silico
            )
        )

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
    return round((n_confirmed/n_silico)*100, 2)

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
    print(f"in silico fragments = {insilico_fragments}")
    print(f"confirmed fragments = {confirmed_fragments}")
    
    # work out percentage found fragments from core series 
    percentage_found_fragments = get_percent_found_fragments(
        insilico_fragments=insilico_fragments,
        confirmed_fragments=confirmed_fragments,
        core_fragment_series=core_fragment_series,
        optional_fragments=optional_fragments,
        exclude_fragments=exclude_fragments
    )

    # initiate list to store missing essential fragments
    missing_essentials = []

    # check that all essential fragments are in in silico fragments
    if essential_fragments:
        essential_fragments = list(
            filter(
                lambda frag: frag in insilico_fragments, essential_fragments
            )
        )
    
    # if any essential fragments are missing, restrict confidence assignment
    # to an upper limit of essential_fragment_cap
        missing_essentials = [
            frag for frag in essential_fragments
            if frag not in confirmed_fragments
        ]

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
            silico_fragments=insilico_fragments,
            confirmed_fragments=confirmed_fragments,
            core_fragment_series=core_fragment_series,
            optional_fragments=optional_fragments,
            excluded_fragments=exclude_fragments
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

def get_Rt_I_binning(
    target_EIC,
    target_Rt,
    Rt_bin
):
    """
    Takes an EIC, a target retention time and target bin for screening
    retention time, and returns the most intense Rt_I pair within that bin
    region

    Args:
        target_EIC (list of lists): list of Rt_I sublists
        target_Rt (float): target retention time
        Rt_bin (float): bin tolerance for returning most intense Rt_I pair -
            most intense data point within range
            (target_Rt-Rt_bin, target_Rt+Rt_bin) will be returned

    Returns:
        Rt_I: Rt_I pair for most intense data point within Rt_bin
    """
    # remove data points that do not fit within Rt_bin
    working_EIC = list(
        filter(
            lambda Rt_I: abs(Rt_I[0]-Rt_bin) >= target_Rt,
            target_EIC
        )
    )

    # if any data is left, return most intense data point
    if working_EIC:
        working_EIC = sorted(
            working_EIC, key=lambda Rt_I: Rt_I[1]
        )
        return working_EIC[-1]

    # return nothing if no data survives Rt_bin filter
    return []

def get_Rt_I_from_ms2_EIC(
    MS1_EIC,
    MS2_EIC,
    ms2_Rt_bin,
    flexible_ms2_rt=True
):
    """
    Takes an MS1 and MS2 EIC for target, and returns Rt_I intensity from MS1
    by matching MS1 EIC to most intense data point in MS2 EIC

    Args:
        MS1_EIC (list of lists): list of [Rt, I] sublists for MS1 EIC
        MS2_EIC (list of lists): list of [Rt, I] sublists for MS2 EIC
        ms2_Rt_bin (float): specifies retention time tolerance for matching
            MS2 retention time to MS2 retention time
        flexible_ms2_rt (bool, optional): specifies whether to widen ms2_Rt_bin
            until matching MS1 EIC [Rt, I] is found. Defaults to True.

    Raises:
        Exception: raise Exception if MS1 EIC is empty

    Returns:
        ms1_Rt_I: MS1 retention time and intensity that most closely matches
            MS2 data within the constraints set
    """
    if len(MS1_EIC) == 0:
        raise Exception('you have an empty MS1 EIC. This is a very serious problem. It\'s bad, and you should feel bad!')

    # get MS2 EIC retention time
    ms2_Rt = sorted(
        MS2_EIC, key=lambda Rt_I: Rt_I[1]
    )[-1][1]

    # set up working copy of MS1 EIC
    working_ms1_EIC = MS1_EIC

    # initate empty list to store assigned MS1 Rt_I
    ms1_Rt_I = []

    # if no flexibility is NOT permitted, return whatever Rt_I data point is
    # found for one round of screening MS1 EIC
    if not flexible_ms2_rt:
        ms1_Rt_I = get_Rt_I_binning(
            target_EIC=working_ms1_EIC,
            target_Rt=ms2_Rt,
            Rt_bin=ms2_Rt_bin
        )

        return ms1_Rt_I

    # if flexibility is permittd, keep doubling Rt_bin for matching MS1 and MS2
    # EICs until an Rt_I match is found and returned
    while not ms1_Rt_I:

        ms1_Rt_I = get_Rt_I_binning(
            target_EIC=working_ms1_EIC,
            target_Rt=ms2_Rt,
            Rt_bin=ms2_Rt_bin
        )

        ms2_Rt_bin = ms2_Rt_bin*2

    return ms1_Rt_I

def return_peaks_EIC(
    EIC,
    n_peaks,
    Rt_bin,
    min_relative_intensity=None
):
    """
    Takes an EIC and returns peaks that fit within retention time limit
    (Rt_bin), starting from most intense peak

    Args:
        EIC (list of lists): list of [Rt, I] sublists
        n_peaks (int): number of peaks to return
        Rt_bin (float): minimum gap between peaks (in retention time)
        min_relative_intensity (float, optional): minimum intensity of peaks
            as a DECIMCAL FRACTION OF MOST INTENSE PEAK. Defaults to None.

    Returns:
        peaks (list of lists): list of [Rt, I] sublists for matching peaks
    """

    # initiate list to store found peaks
    peaks = []

    # create working copy of input EIC and sort by intensity
    working_ms1_EIC = sorted(
        EIC, key=lambda Rt_I: Rt_I[1]
    )

    # set minimum relative intensity to 0 if no value has been specified
    if not min_relative_intensity:
        min_relative_intensity = 0

    # iterate through EIC for specified number of targets
    for i in range(0, n_peaks):

        # find most intense peak in EIC, add to peaks
        peak = working_ms1_EIC[-1]
        peaks.append(peak)

        # return found peaks if all targets are accounted for
        if len(peaks) == n_peaks:
            return peaks

        # remove peaks that are too close to previously found target peaks
        working_ms1_EIC = list(
            filter(
                lambda Rt_I:
                abs(Rt_I[0]-peak[0]) >= Rt_bin
                and Rt_I[1] >= min_relative_intensity*peak[-1],
                working_ms1_EIC
            )
        )

        # ensure EIC remains sorted by intensity
        working_ms1_EIC = sorted(
            working_ms1_EIC, key=lambda Rt_I: Rt_I[1]
        )

    # return list of found targets
    return peaks

def get_Rt_I_from_ms1_EIC(
    EIC,
    Rt_bin,
    backup_Rt_bin=None,
    n_targets=1,
    min_relative_intensity=None
):
    """
    Takes an MS1 EIC and returns matching MS1 peaks that match target(s)
    within constraints set
    Args:
        EIC (list of lists): list of [Rt, I] sublists
        Rt_bin (float): minimum gap between peaks (in retention time)
        backup_Rt_bin (float, optional): backup to Rt_bin: minimum gap between
            peaks if Rt_bin is too stringent for data. Defaults to None.
        n_targets (int, optional): specifies number of target peaks to
            return. Defaults to 1.
        min_relative_intensity (float, optional): minimum relative intensity
            for peaks to be returned expressed as a DECMICAL FRACTION OF MOST
            INTENSE PEAK IN THE EIC. Defaults to None.

    Returns:
        found_targets (list of lists): list of [Rt, I] sublists for found peaks
    """
    # create working copy of EIC and sort by intensity
    working_ms1_EIC = sorted(
        EIC, key=lambda Rt_I: Rt_I[1]
    )

    # return found targets within standard Rt_bin tolerance
    found_targets = return_peaks_EIC(
        EIC=working_ms1_EIC,
        n_peaks=n_targets,
        Rt_bin=Rt_bin,
        min_relative_intensity=min_relative_intensity
    )

    # if back-up Rt_bin has been supplied, use this if required
    if len(found_targets) < n_targets:
        if backup_Rt_bin:
            found_targets = return_peaks_EIC(
                EIC=working_ms1_EIC,
                n_peaks=n_targets,
                Rt_bin=backup_Rt_bin,
                min_relative_intensity=min_relative_intensity
            )

    # return [Rt, I] sublist of peaks found
    return found_targets

def create_EIC_plots(EICs_json):
    """
    This function creates an EIC plot (retention time vs intensity) for
    each sequence in the EICS.json and saves each figure as a .png files.

    Arguments:
        EICs_json (str) -- ".._EICs.json" created by previous post processing
            functions.

    """

    # read EICs.json
    with open(EICs_json, 'r') as EICs_json_data:
        EICs_dict = json.load(EICs_json_data)

    # create pandas dataframe from json data
    EICs_normalized_data = pd.DataFrame.from_dict(
        json_normalize(EICs_dict),
        orient='columns'
    )

    # change dataframe layout (sequence becomes indexed variable)
    EICs_data = pd.melt(EICs_normalized_data)

    # count number of rows
    nrows = len(EICs_data)

    # create new columns for retention times and intensities
    EICs_data["retention times"] = EICs_data["value"]
    EICs_data["intensities"] = EICs_data["value"]

    # separate retention times and intensities to different rows
    for i in range(nrows):
        rt = list(zip(*EICs_data["value"][i]))[0]
        EICs_data["retention times"][i] = rt
        intensity = list(zip(*EICs_data["value"][i]))[1]
        EICs_data["intensities"][i] = intensity

    # remove 'value' column, rename 'variable' to sequence
    EICs_data = EICs_data.drop(columns = 'value')
    EICs_data.rename(columns = {'variable':'sequence'}, inplace=True)

    # create EIC plot for each sequence (retention time vs intensity)
    sns.set_style("ticks")
    sns.set_palette("Set2")
    for i in range(nrows):
        plt.close()
        plt.xlabel("Retention Time")
        plt.ylabel("Intensity")
        plt.title(EICs_data["sequence"][i])
        sns.lineplot(
            (
                EICs_data["retention times"][i]),
            (
                EICs_data["intensities"][i]), data=EICs_data
            )

        # save figure as .png
        plt.savefig(f"{i}_{EICs_data['sequence'][i]}.png")
