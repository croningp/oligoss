import os
import logging
import mzmlripper.extractor as ripper

from bisect import bisect_left, bisect_right

from .extractor_classes import RipperDict
from ..utils.file_io import open_json, write_to_json, return_jsons
from ..utils.run_utils import exception_handler

def mzml_to_json(input_folder, extractor_parameters):
    """ This function takes all mzml files in the input folder and converts them
    to json format using Graham Keenan's ripper if required.

    Args:
        input_folder (str): filepath location of input mzML's
        extractor_parameters (object): full extractor parameters from parameter
        handlers
    """
    if return_jsons(input_folder):
        logging.info('mzml to json conversion not required')
        return return_jsons(input_folder)

    # get list of mzml files in the input folder
    mzml_files = [
        file for file in os.listdir(input_folder) if file.endswith('.mzML')]

    logging.info(f'converting {len(mzml_files)} mzml files to json')

    # convert them to .json format using Graham's ripper
    for mzml in mzml_files:
        ripper.process_mzml_file(
            filename=(os.path.join(input_folder, mzml)),
            out_dir=input_folder,
            rt_units=extractor_parameters.rt_units)

    logging.info('mzml to json conversion complete')

def apply_prefilters(spectra, min_rt, max_rt, min_max_int, min_total_int):
    """ This function filters spectra by retetention time then minimum maximum
    spectra intensity.

    Args:
        spectra (dict): dictionary of format {spectrum_id : {m/z: intensity}}
        min_rt (float): minimum retention time
        max_rt (float): maximum retention time
        min_max_int (float): minimum intensity for the most intense peak
            in the spectrum
        min_total_int (float): minimum intensity for the total intensity
            of the spectrum

    Returns:
        rt_int_filtered_spectra (dict): spectra that has passed the retention
        time and intensity filters
    """
    # filter by retention time
    rt_filtered_spectra = rt_filter(spectra, min_rt, max_rt)

    if not rt_filtered_spectra:
        logging.info('no spectra passed retention time filters')

    # filter by intensity
    rt_int_filtered_spectra = intensity_filter(
        spectra=rt_filtered_spectra,
        min_max_intensity=min_max_int,
        min_total_intensity=min_total_int)

    return rt_int_filtered_spectra

def apply_prefilters_ms2(
        spectra, min_rt, max_rt, min_max_intensity, min_total_intensity,
        ms2_precursors, error, error_units):
    """ This function filters MS2 spectra in the following order: retention time,
        minimum maximum intensity, precursor.

    Args:
        spectra (Dict[str, dict]): dictionary of format
            {spectrum_id : {m/z: intensity}}
        min_rt (float): minimum retention time
        max_rt (float): maximum retention time
        min_max_intensity (float): minimum intensity for the most intense peak
        in the spectrum
        min_total_intensity (float): minimum total intensity for the sum of all
            peaks in the spectrum
        ms2_precursors (List[float]): list of precursor m/z values
        error (float): acceptable difference between target mass and found
            parent mass
        error_units (str): 'ppm' or 'abs'

    Returns:
        filtered_spectra (dict): spectra that has passed the retention time,
        intensity filters and parent mass filters
    """
    # filter by retention time and intensity
    filtered_spectra = apply_prefilters(
        spectra=spectra,
        min_rt=min_rt,
        max_rt=max_rt,
        min_max_int=min_max_intensity,
        min_total_int=min_total_intensity)

    # if ms2, filter by precursors
    filtered_spectra = find_precursors(
        spectra=filtered_spectra,
        ms2_precursors=ms2_precursors,
        error=error,
        error_units=error_units)

    return filtered_spectra

def rt_filter(spectra, min_rt, max_rt):
    """
    This function filters spectra based on retention time. Spectra within the
    minimum and maximum retention time range will be returned.

    Args:
        spectra (Dicr[str, dict]): dictionary of format:
            {spectrum_id : {m/z: intensity}}
        min_rt (float): minimum retention time
        max_rt (float): maximum retention time

    Returns:
        dict: spectra that has passed the retention time filter.
    """
    if not min_rt and not max_rt:
        return spectra

    if not min_rt:
        min_rt = 0

    if not max_rt:
        max_rt = 1E6

    # make ordered lists of retention times and spectra ids
    rt_list = []
    spectrum_id_list = []

    for spectrum_id, spectrum_dict in spectra.items():
        rt_list.append(float(spectrum_dict["retention_time"]))
        spectrum_id_list.append(spectrum_id)

    # find rt values closest to minimum and maximum
    lowest_rt = bisect_right(rt_list, min_rt)
    highest_rt = bisect_left(rt_list, max_rt)

    # get the list of spectra within these limits, build the filtered dictionary
    spectra_list = spectrum_id_list[lowest_rt:highest_rt]

    return {key: spectra[key] for key in spectra_list}

def intensity_filter(spectra, min_max_intensity, min_total_intensity):
    """ This function filters based on the minimum maximum intensity and the
    minimum total intensity of the spectra.
    Spectra which pass intensity intensity thresholds will be returned.

    Args:
        spectra (dict): dictionary of format {spectrum_id : {m/z: intensity}}
        min_max_intensity (float): minimum intensity for the most intense peak
        in the spectrum
        min_total_intensity (float): minimum intensity for the total intensity
        of the spectrum

    Returns:
        total_intensity_filtered (dict): spectra that has passed the intensity
        filter.
    """
    # filter by minimum maximum intensity
    max_intensity_filtered = min_max_intensity_filter(
        spectra=spectra,
        min_max_intensity=min_max_intensity)

    # filter by minimum total intensity
    total_intensity_filtered = min_total_intensity_filter(
        spectra=max_intensity_filtered,
        min_total_intensity=min_total_intensity)

    return total_intensity_filtered


def min_max_intensity_filter(spectra, min_max_intensity):
    """ This function filters based on the maximum intensity in the spectra.
    Spectra which maximum intensity exceeds the limit will be returned.

    Args:
        spectra (dict): dictionary of format {spectrum_id : {m/z: intensity}}
        min_max_intensity (float): minimum intensity for the most intense peak
            in the spectrum

    Returns:
        max_int_filtered (dict): spectra that has passed the intensity filter.
    """
    # return unfiltered spectra if no filter spectified
    if min_max_intensity is None:
        return spectra

    max_int_filtered = {}

    # find whether maximum intensity in each spectrum meets minimum limit
    for spectrum_id, spectrum_dict in spectra.items():
        max_intensity = max(
            [float(x) for x in list(spectrum_dict.values()) if type(x) == int])

        if float(max_intensity) >= float(min_max_intensity):
            max_int_filtered[spectrum_id] = spectrum_dict

    if not max_int_filtered:
        logging.info('no spectra passed maximum intensity filters')

    return max_int_filtered

def min_total_intensity_filter(spectra, min_total_intensity=None):
    """ This function filters based on the total intensity in the spectra.
    Spectra which total intensity exceeds the limit will be returned.

    Args:
        spectra (dict): dictionary of format {spectrum_id : {m/z: intensity}}
        min_total_intensity (float optional): minimum intensity for the total
        intensity of the spectrum. Defaults to None.

    Returns:
        total_int_filtered: spectra that has passed the intensity filter
    """
    # return unfiltered spectra if no filter spectified
    if not min_total_intensity:
        return spectra

    total_int_filtered = {}

    # find whether maximum intensity in each spectrum meets minimum limit
    for spectrum_id, spectrum_dict in spectra.items():
        total_intensity = sum(
            [float(x) for x in list(spectrum_dict.values()) if type(x) == int])

        if float(total_intensity) >= float(min_total_intensity):
            total_int_filtered[spectrum_id] = spectrum_dict

    if not total_int_filtered:
        logging.info('no spectra passed total intensity filters')

    return total_int_filtered

def find_precursor(
    ms2_precursor,
    error,
    error_units,
    connection=None,
    ripper_name=None,
    spectra=None
):
    """ This function returns MS2 spectra which parent ion matches the target
    mass.
    Args:
        ms2_precursor (float): precursor m/z value
        error (float): acceptable difference between target mass and found
            parent mass
        error_units (str): 'abs' or 'ppm'.
        connection (str): polymersoup localhost connection string.
            Defaults to None as not needed during prefiltering.
        ripper_name (str): unique str tag for original ripper file.
            Defaults to None as not needed during prefiltering.
        spectra (dict]): dictionary of multiple spectra, format:
            {spectrum_id : {m/z: intensity}}. Defaults to None.
            Spectra needed for prefiltering but not for later database queries.

    Returns:
        precursor_filtered dict: dict of MS2 spectra whos parent ion matches the
        target mass
    """
    precursor_filtered = {}

    # make sure error is abs
    if error_units == 'ppm':
        error = (float(ms2_precursor) / 1E6 * error)

    # get list of target precursor masses
    target = [float(ms2_precursor) - error, float(ms2_precursor) + error]

    # if spectra is provided, including empty spectra, ie for prefiltering
    # functions iterate through & return filtered spectra with precursor matches
    if spectra or spectra == {}:

        # find whether precursor ion mass matches a target
        for spectrum_id, spectrum_dict in spectra.items():
            parent_mass = float(spectrum_dict["parent"])

            # iterate through targets checking if the parent mass is a match
            if parent_mass >= target[0] and parent_mass <= target[1]:
                precursor_filtered[spectrum_id] = spectrum_dict

    # if no spectra is provided, ie for ms2 screening functions,
    # find spectra matches in database
    else:

        # return spectra with parent masses are within the target range
        # $gt = greater than, $lt = less than for range query
        ms2_spectra = list(
            connection.polymersoup[f"{ripper_name}_ms2_data"].find(
                {"parent": {"$gt": target[0], "$lt": target[1]}}))

        # add all 'spectrum' parts of each match to a dict
        precursor_filtered = {
            entry["spectrum_id"]: entry["spectrum"] for entry in ms2_spectra}

        connection.close()

    return precursor_filtered

def find_precursors(
    ms2_precursors,
    error,
    error_units,
    spectra=None,
    connection=None,
    ripper_name=None
):
    """ This function returns MS2 spectra which parent ion matches the target
    mass for all targets.

    Args:
        ms2_precursors (List[float]): list of MS2 precrursor m/z values
        error (float): acceptable difference between target mass and found
            parent mass. Can either be in absolute units (u) or relative (ppm).
        error_units (str): 'ppm' or 'abs'.
        connection (str): polymersoup localhost connection string
        ripper_name (str): unique str tag for original ripper file.
        spectra (dict]): dictionary of multiple spectra, format:
            {spectrum_id : {m/z: intensity}}.
            Defaults to None as not needed during ms2 screening when spectra is
            retrieved from database queries.

    Returns:
        precursor_filtered (dict): dict of MS2 spectra whos parent ion matches
            any target mass
    """
    precursor_filtered = {}

    # if ms2_precursor argument isn't given for prefiltering, return spectra
    if not ms2_precursors and spectra:
        return spectra

    for ms2_precursor in ms2_precursors:
        precursor_filtered.update(find_precursor(
            spectra=spectra, ms2_precursor=ms2_precursor, error=error,
            error_units=error_units, connection=connection,
            ripper_name=ripper_name))

    return precursor_filtered

def min_ms2_peak_abundance_filter(
        spectra, peak_list, error, min_ms2_peak_abundance=None):
    """ This function checks all spectra meet the minimum ms2 peak abundance.
    Each spectrum is screened for the most intense peak from the peak list that
    is associated with the sequence. The intensity of this peak is compared with
    the most intense peak in the spectrum (not-related to sequence) and the
    spectrum is returned if it exceeds the min ms2 peak abundance threshold.

    Args:
        spectra (Dict[str, dict]): dictionary of MS2 spectra to be filtered.
        peak_list (List[float]): list of all peaks (floats) associated with a
        given sequence.
        error (float): absolute error tolerance for matching target ions to
        masses found in spectra
        min_ms2_peak_abundance (float, optional): minimum percentage of ms2
        peak abundance. Defaults to None.

    Returns:
        min_ms2_peak_filtered (Dict[str, dict]): dictionary of spectra that
            have passed the filter.
    """
    # return unfiltered spectra if no filter spectified
    if not min_ms2_peak_abundance:
        return spectra

    min_ms2_peak_filtered = {}

    # find whether peak associated with sequence passes abundance limit
    for spectrum_id, spectrum_dict in spectra.items():

        # find most intense peak in spectrum
        max_intensity = max(
            [float(x) for x in list(spectrum_dict.values()) if type(x) == int])

        # find most intense peak in spectrum associated to sequence
        highest_int_frag = 0
        for mass in peak_list:
            matches = match_mass(
                spectrum=spectrum_dict,
                mass_range=[mass - error, mass + error])
            for match in matches:
                if spectrum_dict[str(match)] > highest_int_frag:
                    highest_int_frag = spectrum_dict[str(match)]

        ms2_peak_abundance = (
            float(highest_int_frag) / float(max_intensity)) * 100

        # if ms2 peak abundance is higher than threshold, pass the spectrum
        if ms2_peak_abundance >= min_ms2_peak_abundance:
            min_ms2_peak_filtered[spectrum_id] = spectrum_dict

    return min_ms2_peak_filtered

def match_mass(spectrum, mass_range):
    """ This function takes a mass range for a single target and finds matches
    within the spectrum.

    Args:
        spectrum (dict): single spectrum, must contain 'mass_list'.
        mass_range (list): format: [minimum mass, maximum mass]

    Returns:
        matches (list[float]): list of mass matches from spectrum that fall
        within limit.
    """
    matches = []

    # find closest values to min and max mass in mass list
    mass_list = spectrum['mass_list']
    min_match = bisect_left(mass_list, mass_range[0])
    max_match = bisect_right(mass_list, mass_range[1])

    # find matches within mass range
    try:
        mass_matches = [mass_list[min_match]]

    # if there's no matches (index exceeds the list), return empty list
    except IndexError:
        return matches

    if min_match != max_match:
        mass_matches = mass_list[min_match:max_match]

    if mass_matches:
        matches = [
            format(m, '.4f') for m in mass_matches if (
                m >= mass_range[0]) and (m <= mass_range[1])]

    return matches

@exception_handler(fatal=True)
def prefilter(ripper_file, extractor_parameters):
    """ This function applies all prefilters (retention time, minimum maximum
    intensity threshold, minimum total intensity threshold, precursor mass) to
    all spectra for a single ripper file.
    The ripper json is converted to an object and its spectra are filtered by
    the above parameters, starting with MS1 then moving onto MS2.

    Args:
        ripper_file (str): full filepath to ripper json
        extractor_parameters (object): full extractor parameters from parameter
        handlers

    Returns:
        full_filepath (str): full filepath to filtered ripper json
    """
    # check for pre screen filters
    pre_screen_filters = extractor_parameters.pre_screen_filters

    if not pre_screen_filters:
        return ripper_file

    #  define pre-screen filters
    if "min_ms1_max_intensity" in pre_screen_filters:
        min_max_ms1 = pre_screen_filters["min_ms1_max_intensity"]
    else:
        min_max_ms1 = 0

    if "min_ms2_max_intensity" in pre_screen_filters:
        min_max_ms2 = pre_screen_filters["min_ms2_max_intensity"]
    else:
        min_max_ms2 = 0

    if "min_ms1_total_intensity" in pre_screen_filters:
        min_total_ms1 = pre_screen_filters["min_ms1_total_intensity"]
    else:
        min_total_ms1 = 0

    if "min_ms2_total_intensity" in pre_screen_filters:
        min_total_ms2 = pre_screen_filters["min_ms2_total_intensity"]
    else:
        min_total_ms2 = 0

    if "min_rt" in pre_screen_filters:
        min_rt = pre_screen_filters["min_rt"]
    else:
        min_rt = 0

    if "max_rt" in pre_screen_filters:
        max_rt = pre_screen_filters["max_rt"]
    else:
        max_rt = 0

    if "precursors" in pre_screen_filters:
        ms2_precursors = pre_screen_filters["precursors"]
    else:
        ms2_precursors = None

    # open ripper json and convert dict to object
    ripper_data = open_json(ripper_file)
    ripper = RipperDict(ripper_data)

    # apply all pre-filters to each ms1 spectrum
    ripper.spectra["ms1"] = apply_prefilters(
        spectra=ripper.ms1,
        min_rt=min_rt,
        max_rt=max_rt,
        min_max_int=min_max_ms1,
        min_total_int=min_total_ms1)

    # check format of min ms2 max intensity is absolute
    if float(min_max_ms2) <= 1:
        min_max_ms2 = float(min_max_ms1) * float(min_max_ms2)

    # save to json and return filepath
    output_folder = os.path.join(os.path.dirname(ripper_file), "prefiltered")
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    full_filepath = os.path.join(
        output_folder, os.path.basename(ripper_file))

    # apply all pre-filters to each ms2 spectrum
    ripper.spectra["ms2"] = apply_prefilters_ms2(
        spectra=ripper.ms2,
        min_rt=min_rt,
        max_rt=max_rt,
        min_max_intensity=min_max_ms2,
        min_total_intensity=min_total_ms2,
        ms2_precursors=ms2_precursors,
        error=extractor_parameters.error,
        error_units=extractor_parameters.error_units)

    write_to_json(ripper.spectra, full_filepath)

    #  explicitly reassign ripper to NoneType to avoid memory leaks. NOTE: I am
    #  aware this looks daft - if you want to discuss it let me know. David.
    ripper = {}

    return full_filepath

@exception_handler(fatal=True, verbose=True)
def prefilter_all(rippers, extractor_parameters):
    """ This function applies all prefilters (retention time, minimum maximum
    intensity threshold, minimum total intensity, precursor mass) to all ripper
    file spectra.

    Args:
        rippers (str): full filepath to ripper json
        extractor_parameters (object): full extractor parameters from parameter
        handlers

    Returns:
        rippers (List[str]): list of full filepaths to all filtered rippers
    """

    if extractor_parameters.pre_screen_filters:
        logging.info(
            f"applying pre-screen filters to {len(rippers)} ripper JSONs")

        filtered_rippers = []

        for i, ripper_file in enumerate(rippers):
            logging.info(
                f'applying prefilters to {ripper_file} (file {i + 1} of\n'
                f'{len(rippers)})'.rstrip('\n'))

            # pre filter spectra
            filtered_rippers.append(prefilter(
                ripper_file=ripper_file,
                extractor_parameters=extractor_parameters
            ))

        return filtered_rippers

    logging.info('prefilters not specified')
    return rippers
