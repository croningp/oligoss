"""
This file contains functions for extracting sequence-specific data from
spectra

"""
from .extractor_helpers import *

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
        precursors ([type]): [description]
        fragment_dict ([type]): [description]
        peak_list ([type]): [description]
        ripper_dict ([type]): [description]
        err ([type]): [description]
        err_abs (bool, optional): [description]. Defaults to True.
        ms_level (int, optional): [description]. Defaults to 2.
        min_annotated_peak_assignment (int, optional): [description]. Defaults to 90.
        essential_signatures (list, optional): [description]. Defaults to [].
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

    # iterate through spectra and find most intense peak
    for spectrum in ripper_dict:

        most_intense_peak = find_most_intense_peak_spectrum



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
