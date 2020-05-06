import os
import mzmlripper.extractor as ripper
from bisect import bisect_left, bisect_right
from extractor_classes import ripper_dict
from general_functions import logging, open_json, write_to_json, return_jsons

def mzml_to_json(input_folder, extractor_parameters):
    """
    This function takes all mzml files in the input folder and converts them
    to json format using Graham Keenan's ripper if required.

    Arguments:
        input_folder {str} -- filepath location of input mzML's
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
            out_dir=input_folder)

    return logging.info('mzml to json conversion complete')

def apply_prefilters(spectra, min_rt, max_rt, min_max_int, min_total_int):
    """
    This function filters spectra by retetention time then minimum maximum
    spectra intensity.

    Arguments:
        spectra {dict} -- dictionary of format {spectrum_id : {m/z: intensity}}
        min_rt {float} -- minimum retention time
        max_rt {float} -- maximum retention time
        min_max_intensity {float} -- minimum intensity for the most intense peak
            in the spectrum
        min_total_intensity {float} -- minimum intensity for the total intensity
            of the spectrum
        ms_level {int} -- MS1 or MS2

    Returns:
        dict -- spectra that has passed the retention time and intensity filters
    """
    # filter by retention time
    rt_filtered_spectra = rt_filter(spectra, min_rt, max_rt)

    if not rt_filtered_spectra:
        logging.info(f'no spectra passed retention time filters')

    # filter by intensity
    rt_int_filtered_spectra = intensity_filter(
        spectra=rt_filtered_spectra,
        min_max_intensity=min_max_int,
        min_total_intensity=min_total_int)

    return rt_int_filtered_spectra

def apply_prefilters_ms2(
        spectra, min_rt, max_rt, min_max_intensity, min_total_intensity,
        ms2_precursors, error):
    """
    This function filters MS2 spectra in the following order: retention time,
        minimum maximum intensity, precursor.

    Arguments:
        spectra {dict} -- dictionary of format {spectrum_id : {m/z: intensity}}
        min_rt {float} -- minimum retention time
        max_rt {float} -- maximum retention time
        min_max_intensity {float} -- minimum intensity for the most intense peak
            in the spectrum
        ms_level {int} -- MS1 or MS2 (1 or 2)
        ms2_precursors {list} -- list or nested list of format:
            [precursor string, precursor mass]
        error {float} -- acceptable difference between target mass and found
            parent mass

    Returns:
        dict -- spectra that has passed the retention time, intensity filters
            and parent mass filters
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
        error=error)

    return filtered_spectra

def rt_filter(spectra, min_rt, max_rt):
    """
    This function filters spectra based on retention time. Spectra within the
    minimum and maximum retention time range will be returned.

    Arguments:
        spectra {dict} -- dictionary of format {spectrum_id : {m/z: intensity}}
        min_rt {float} -- minimum retention time
        max_rt {float} -- maximum retention time

    Returns:
        dict -- spectra that has passed the retention time filter.
    """

    if min_rt is None or max_rt is None:
        return spectra

    if min_rt is None:
        min_rt = 0

    if max_rt is None:
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
    """
    This function filters based on the minimum maximum intensity and the
    minimum total intensity of the spectra.
    Spectra which pass intensity intensity thresholds will be returned.

    Arguments:
        spectra {dict} -- dictionary of format {spectrum_id : {m/z: intensity}}
        min_max_intensity {float} -- minimum intensity for the most intense peak
            in the spectrum
        min_total_intensity {float} -- minimum intensity for the total intensity
            of the spectrum

    Returns:
        dict -- spectra that has passed the intensity filter.
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
    """
    This function filters based on the maximum intensity in the spectra.
    Spectra which maximum intensity exceeds the limit will be returned.

    Arguments:
        spectra {dict} -- dictionary of format {spectrum_id : {m/z: intensity}}
        min_max_intensity {float} -- minimum intensity for the most intense peak
            in the spectrum

    Returns:
        dict -- spectra that has passed the intensity filter.
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
        logging.info(f'no spectra passed maximum intensity filters')

    return max_int_filtered

def min_total_intensity_filter(spectra, min_total_intensity=None):
    """
    This function filters based on the total intensity in the spectra.
    Spectra which total intensity exceeds the limit will be returned.

    Arguments:
        spectra {dict} -- dictionary of format {spectrum_id : {m/z: intensity}}
        min_total_intensity {float} -- minimum intensity for the total intensity
            of the spectrum

    Returns:
        dict -- spectra that has passed the intensity filter.
    """
    # return unfiltered spectra if no filter spectified
    if min_total_intensity is None:
        return spectra

    total_int_filtered = {}

    # find whether maximum intensity in each spectrum meets minimum limit
    for spectrum_id, spectrum_dict in spectra.items():
        total_intensity = sum(
            [float(x) for x in list(spectrum_dict.values()) if type(x) == int])

        if float(total_intensity) >= float(min_total_intensity):
            total_int_filtered[spectrum_id] = spectrum_dict

    if not total_int_filtered:
        logging.info(f'no spectra passed total intensity filters')

    return total_int_filtered

def find_precursor(spectra, ms2_precursor, error):
    """
    This function returns MS2 spectra which parent ion matches the target mass.

    Arguments:
        spectra {dict} -- dictionary of format {spectrum_id : {m/z: intensity}}
        ms2_precursor {list} -- [ms2 precursor string, ms2 precursor mass]
        error {float} -- acceptable difference between target mass and found
            parent mass

    Returns:
        dict -- dict of MS2 spectra whos parent ion matches the target mass
    """

    precursor_filtered = {}

    # get list of target precursor masses
    target = [float(ms2_precursor[1]) - error, float(ms2_precursor[1]) + error]

    # find whether precursor ion mass matches a target
    for spectrum_id, spectrum_dict in spectra.items():
        parent_mass = float(spectrum_dict["parent"])

        # iterate through targets checking if the parent mass is a match
        if parent_mass >= target[0] and parent_mass <= target[1]:
            precursor_filtered[spectrum_id] = spectrum_dict

    return precursor_filtered

def find_precursors(spectra, ms2_precursors, error):
    """
    This function returns MS2 spectra which parent ion matches the target mass
    for all targets.

    Arguments:
        spectra {dict} -- dictionary of format {spectrum_id : {m/z: intensity}}
        ms2_precursor {list or nested list} -- format:
            [ms2 precursor string, ms2 precursor mass]
        error {float} -- acceptable difference between target mass and found
            parent mass

    Returns:
        dict -- dict of MS2 spectra whos parent ion matches any target mass
    """

    if ms2_precursors is None:
        return spectra

    # check if there are multiple ms2 precursors (nested list)
    if type(ms2_precursors[0]) != list:
        return find_precursor(
            spectra=spectra,
            ms2_precursor=ms2_precursors,
            error=error)

    precursor_filtered = {}
    for ms2_precursor in ms2_precursors:
        precursor_filtered = precursor_filtered.update(find_precursor(
            spectra=spectra, ms2_precursor=ms2_precursor, error=error))

    return precursor_filtered

def prefilter(ripper_file, extractor_parameters, input_folder):
    """
    This function applies all prefilters (retention time, minimum maximum
    intensity threshold, minimum total intensity threshold, precursor mass) to
    all spectra for a single ripper file.
    The ripper json is converted to an object and it's spectra are filtered by
    the above parameters, starting with MS1 then moving onto MS2.

    Arguments:
        ripper_file {str} -- full filepath to ripper json
        extractor_parameters {object} -- full extractor parameters
        input_folder {str} -- full filepath to ripper json input folder

    Returns:
        str -- full filepath to filtered ripper json
    """
    # open ripper json and convert dict to object
    ripper_data = open_json(ripper_file)
    ripper = ripper_dict(ripper_data)
    min_max_ms1 = extractor_parameters.filter_min_ms1_max_intensity
    min_max_ms2 = extractor_parameters.filter_min_ms2_max_intensity
    min_total_ms1 = extractor_parameters.filter_min_ms1_total_intensity
    min_total_ms2 = extractor_parameters.filter_min_ms2_total_intensity

    # make sure error is absolute
    if extractor_parameters.error_units == 'ppm':
        error = float(extractor_parameters.error) / 1E6

    else:
        error = float(extractor_parameters.error)

    # apply all pre-filters to each ms1 spectrum
    ripper.spectra["ms1"] = apply_prefilters(
        spectra=ripper.ms1,
        min_rt=extractor_parameters.filter_min_rt,
        max_rt=extractor_parameters.filter_max_rt,
        min_max_int=min_max_ms1,
        min_total_int=min_total_ms1)

    # check format of min ms2 max intensity is absolute
    if float(min_max_ms2) <= 1:
        min_max_ms2 = float(min_max_ms1) * float(min_max_ms2)

    # apply all pre-filters to each ms2 spectrum
    ripper.spectra["ms2"] = apply_prefilters_ms2(
        spectra=ripper.ms2,
        min_rt=extractor_parameters.filter_min_rt,
        max_rt=extractor_parameters.filter_max_rt,
        min_max_intensity=min_max_ms2,
        min_total_intensity=min_total_ms2,
        ms2_precursors=extractor_parameters.precursors,
        error=error)

    # save to json and return filepath
    output_folder = os.path.join(input_folder, "prefiltered")
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    ripper_filename = ripper_file.split('\\')[-1]
    full_filepath = os.path.join(
        output_folder, f"prefiltered_{ripper_filename}")
    write_to_json(ripper.spectra, full_filepath)

    return full_filepath

def prefilter_all(rippers, extractor_parameters, input_folder):
    """ This function applies all prefilters (retention time, minimum maximum
    intensity threshold, minimum total intensity, precursor mass) to all ripper
    file spectra.

    Arguments:
        ripper_file {str} -- full filepath to ripper json
        extractor_parameters {object} -- full extractor parameters
        input_folder {str} -- full filepath to ripper json input folder

    Returns:
        list -- list of full filepaths to all filtered rippers
    """

    if extractor_parameters.pre_screen_filter:

        filtered_rippers = []

        for ripper_file in rippers:
            logging.info(f'applying prefilters to {ripper_file}')

            # pre filter spectra
            filtered_rippers.append(str(prefilter(
                ripper_file=ripper_file,
                extractor_parameters=extractor_parameters,
                input_folder=input_folder)))

        return filtered_rippers
    logging.info(f'prefilters not specified')
    return rippers
