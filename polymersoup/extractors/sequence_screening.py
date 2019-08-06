"""
This file contains functions for extracting sequence-specific data from
spectra

"""
from .run_extractors import *

def generate_EIC(
    ions,
    ms_level,
    ripper_dict,
    err,
    err_abs=True,
    min_max_intensity=None,
    min_total_intensity=None
):
    """
    Takes a list of ions (m/z values) and generates a combined extracted ion
    chromatogram (EIC) of those ions from mass spectra in mzml ripper format
    (ripper_dict)

    Args:
        ions (list of floats): list of target ion m/z values
        ms_level (int): MS level at which to screen for ions
        ripper_dict (dict): mzml ripper mass spectrometry data dictionary
        err (float): error tolerance for matching target ions to masses found
                in spectra; units either in absolute mass units or ppm
        err_abs (bool, optional): specifies whether err is in absolute mass
                units or ppm; if True, err is in absolute mass units.
                Defaults to True.
        min_max_intensity (float, optional): minimum maximum intensity of EIC
                for it to be returned. Defaults to None.
        min_total_intensity (float, optional): maximum total intensity of EIC
                for it to be returned. Defaults to None.

    Returns:
        EIC (list of lists): extracted ion chromatogram in format: [[Rt, I]..]
                where Rt = retention time (float) and I = intensity of all
                species that match target ions at that retention time
    """

    # extract spectra at target ms level
    ripper_dict = ripper_dict[f'ms{ms_level}']

    # initiate list to store EIC
    EIC = []

    # iterate through spectra, and list masses found in each spectrum
    for spectrum in ripper_dict:
        masses = [float(mass) for mass in spectrum['mass_list']]

        # initiate list for storing masses that match target ions
        matches = []

        # iterate through targets, and add matches as they are found for each
        # target to matches list
        for ion in ions:
            matches.extend(find_target(
                ion,
                masses,
                err,
                err_abs)
            )

        # get intensity of all matches in spectrum
        intensity = sum(
            [
                float(spectrum[f'{match}'])
                for match in matches]

        )

        # retrieve spectrum's retention time
        retention_time = float(spectrum['retention_time'])

        # add retention time, intensity as a sublist to EIC
        EIC.append([retention_time, intensity])

    # sort EIC by intensity
    # Rt_I = sublists of [retention_time, intensity] in EIC
    EIC.sort(key = lambda Rt_I: Rt_I[1])

    # check for specified minimum max intensity threshold. If this is specified
    # and most intense point on EIC does not exceed this, return blank EIC
    if min_max_intensity:
        if EIC[-1][1] < min_max_intensity:
            EIC = []

    # check for specified minimum total intensity threshold. If this is
    # specified and sum of EIC intensity does not exceed this, return blank EIC
    if min_total_intensity:
        if sum([Rt_I[1] for Rt_I in EIC]) < min_total_intensity:
            EIC = []

    # return EIC
    return EIC

def generate_MS1_EICs_sequence_dict(
    silico_dict,
    ripper_dict,
    err,
    err_abs,
    min_max_intensity,
    min_total_intensity
):
    """
    Takes a silico_dict of sequences and corresponding MS1 ions, and outputs
    MS1 EICs for sequences

    Args:
        silico_dict ([type]): [description]
        ripper_dict ([type]): [description]
        err ([type]): [description]
        err_abs ([type]): [description]
        min_max_intensity ([type]): [description]
        min_total_intensity ([type]): [description]
    """
    for sequence in silico_dict:
        if type(silico_dict[sequence]) == dict:
            masses = silico_dict[sequence]["MS1"]
        else:
            masses = silico_dict[sequence]

        silico_dict[sequence] = [float(mass) for mass in masses]

        sequence_MS1_EIC = generate_EIC(
            masses,
            1,
            ripper_dict,
            err,
            err_abs,
            min_max_intensity,
            min_total_intensity
        )

def confirm_fragments_sequence(
    precursors,
    fragment_dict,
    peak_list,
    ripper_dict,
    err,
    err_abs=True,
    ms_level=2,
    min_annotated_peak_assignment=90,
    essential_signatures=[]
):
    """
    This function takes a list of precursor ions, fragment dictionary and
    returns list of confirmed fragments associated with precursors

    Args:
        precursors (list of floats): list of m/z values for sequence precursors
        fragment_dict (dict): dictionary of fragments and associated masses
            in insilico fragment dict format
        peak_list (list of floats): list of ALL MS1 and MS2 / MSn masses
            associated with sequence
        ripper_dict (dict): dict of mass spec data in mzml ripper format
        err (float): error tolerance used when matching masses; units can
            either be in absolute mass units or ppm
        err_abs (bool, optional): specifies whether err units are in absolute
            mass units or ppm; if True, absolute mass unit is used.
            Defaults to True.
        ms_level (int, optional): specifies MS level for screening fragments
            and precursors. Defaults to 2.
        min_annotated_peak_assignment (int, optional): minimum % intensity of
            most intense peak in a given spectrum that matches a mass in the
            peak list. Defaults to 90.
        essential_signatures (list, optional): list of signatures that MUST
            be found in a spectrum for fragments to be confirmed from that
            spectrum. Defaults to [].
    Returns:
        confirmed_fragments (list of strings): list of fragment ids for
            fragments in fragment_dict that have been confirmed as present in
            mass spectrum data from ripper_dict
    """
    # remove spectra at fragmentation ms level that do not have matching
    # precursors
    ripper_dict = filter_parent_ion(
        ripper_dict,
        precursors,
        err,
        err_abs,
        ms_level)

    # pull out spectra at fragmentation ms level
    ripper_dict = ripper_dict[f'ms{ms_level}']

    # initiate list to store confirmed fragment ids
    confirmed_fragments = []

    # iterate through spectra
    for spectrum in ripper_dict:

        # calculate annotated_peak_assignemnt from spectrum, peak_list
        annotated_peak_assignment = find_maximum_annotated_peak_intensity(
            spectrum,
            peak_list,
            err,
            err_abs
        )

        # if annotated_peak_assignment is below specified lower threshold,
        # pass over spectrum without searching for fragments
        if annotated_peak_assignment < min_annotated_peak_assignment:
            pass

        # if annotated_peak_assignment is equal to or exceeds threshold, search
        # spectrum for fragments provided in fragment_dict
        else:

            # add fragments found in spectrum to list of confirmed fragments
            # associated with sequence
            confirmed_fragments.extend(
                find_fragments_spectrum(
                    spectrum,
                    fragment_dict,
                    err,
                    err_abs,
                    essential_signatures
                )
            )
    # remove any duplicates from confirmed fragment list
    confirmed_fragments = list(set(confirmed_fragments))

    # finally, return list of confirmed fragments
    return confirmed_fragments

def find_most_intense_peak_spectrum(
    spectrum
):
    """
    Takes a spectrum in mzml_ripper format and returns tuple of most intense
    peak with m/z and associated intensity in format (m/z, I) where m/z =
    m/z value of most intense peak and I = intensity of most intense peak

    Args:
        spectrum (dict): single spectrum dictionary in mzml ripper format

    Returns:
        most_intense (tuple): tuple of peak m/z and intensity in format
                    (m/z, I)
    """
    # get list of mass (m/z), intensity tuples for all peaks in spectrum
    masses = [
        (float(mass), float(spectrum[mass]))
        for mass in spectrum['mass_list']
    ]

    # sort peaks by most intense
    masses.sort(lambda x: x[1])

    # retrieve most intense m/z in format (m/z, I) where m/z = peak m/z and I
    # = intensity of peak with m/z
    most_intense = masses[-1]

    return most_intense

def find_most_intense_matching_peak(
    spectrum,
    peak_list,
    err,
    err_abs):
    """
    takes a list of theoretical peaks associated with a sequence, a spectrum
    containing observed peaks and associated intensities, and returns the
    intensity of the most intense matching peak

    Args:
        spectrum (dict): spectrum dict in mzml ripper dict format
        peak_list (list): list of theoretical m/z values associated with a
                        sequence
        err (float):
        err_abs ([type]): [description]

    Returns:
        top_intensity (float): intensity of most intense matching peak
    """

    # set top matching intensity to 0; this will be reset if matching peaks
    # are found with intensity > 0
    top_intensity = 0

    # iterate through theoretical peaks associated with sequence
    for peak_mass in peak_list:

        # find real peaks in spectrum that match theoretical sequence peak
        # mass
        matches = find_target(
            peak_mass,
            [float(mass) for mass in spectrum['mass_list']],
            err,
            err_abs
        )

        # get intensities of all matching peaks, sorted in ascending order
        intensities = sorted(
            [
                float(spectrum[f'{match}'])
                for match in matches
            ]
        )

        # check whether most intense match is more intense than previous most
        # intense match; if so, reset top_intensity
        if intensities[-1] > top_intensity:
            top_intensity = intensities[-1]

    return top_intensity

def find_maximum_annotated_peak_intensity(
    spectrum,
    peak_list,
    err,
    err_abs
):
    """
    Searches a spectrum in mzml ripper dict format for peaks (m/z values) in
    a peak list and returns % intensity of most intense matching peak relative
    to most intense peak in the spectrum

    Args:
        spectrum (dict): spectrum dict in mzml ripper format
        peak_list (list of floats): list of peaks (m/z values), typically peaks
                    associated with a target sequence
        err (float): error tolerance in matching peak masses - can be in
                    absolute mass units or ppm
        err_abs (bool): specifies whether err units are in absolute mass units
                    or ppm; if True, units = absolute mass units

    Returns:
        float: maximum annotated peak intensity in % of most intense peak
    """

    # find intensity of most intense peak in spectrum
    most_intense_peak = find_most_intense_peak_spectrum(spectrum)[-1]

    # find intensity of most intense spectrum peak that matches something in
    # the peak list
    most_intense_matching_peak = find_most_intense_matching_peak(
        spectrum,
        peak_list,
        err,
        err_abs
    )

    # return maximum annotated peak intensity in % of most intense peak in
    # spectrum
    if most_intense_matching_peak > 0:
        return (most_intense_matching_peak/most_intense_peak) * 100
    else:
        return 0

def find_fragments_spectrum(
    spectrum,
    fragment_dict,
    err,
    err_abs,
    essential_signatures=[]
):
    """
    Takes an MS2 / MSn fragment dict, MS2 / MSn spectrum and returns list of
    fragments with matching ions in the spectrum

    Args:
        spectrum (dict): spectrum dict in mzml ripper format
        fragment_dict (dict): fragment dictionary in standard insilico format
        err (float): error tolerance for matching fragment masses; units can
                    either be absolute mass units or ppm
        err_abs (bool): specifies whether err units are in absolute mass
                    units or ppm; if True, err is in absolute mass units
        essential_signtures (list): list of signature fragment ids that must
                    be found in spectrum if any fragments are to be confirmed
                    from said spectrum (defaults to [])

    Returns:
        confirmed_fragments (list of strings): list of fragment ids from
                    input fragment_dict that have one or more matching masses
                    in the spectrum
    """
    # initiate list to store confirmed fragments
    confirmed_fragments = []

    # retrieve masses from spectrum
    spectrum_masses = [float(mass) for mass in spectrum['mass_list']]

    # iterate through SIGNATURE FRAGMENTS
    for signature, signature_masses in fragment_dict['signatures'].items():

        # retrieve masses in spectrum that match fragment target masses
        matches = find_multiple_targets(
            signature_masses,
            spectrum_masses,
            err,
            err_abs
        )

        # if any spectrum masses are a match for fragment masses, add
        # fragment id to list of confirmed fragments
        if matches:
            confirmed_fragments.append(signature)

    # check whether there are any essential signatures, i.e. signatures that
    # MUST be found in the spectrum for fragments to be confirmed
    if essential_signatures:

        # if any essential signatures are not found, return empty list of
        # confirmed fragments
        for signature in essential_signatures:
            if signature not in confirmed_fragments:
                return []

    # iterate through NON-SIGNATURE fragments
    for fragment, masses in fragment_dict.items():
        if fragment == 'signatures':
            pass

        # retrieve masses in spectrum that match fragment target masses
        matches = find_multiple_targets(
            masses,
            spectrum_masses,
            err,
            err_abs
        )

        # if any spectrum masses are a match for fragment masses, add
        # fragment id to list of confirmed fragments
        if matches:
            confirmed_fragments.append(fragment)

    # return list of fragments found in spectrum
    return confirmed_fragments
