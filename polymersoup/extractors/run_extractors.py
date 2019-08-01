from .extractor_helpers import *
from ..parameter_handlers import *

def pre_filter_rippers(
    ripper_folder='C:/Users/group/polymermassspec/Examples/ripper_jsons',
    parameters_file='C:/Users/group/PolymerMassSpec/Examples/InputParams_Test.json'
):
    """
    Takes a folder containing mzml_ripper jsons, input parameters file for
    data extraction and creates new mzml_ripper_jsons of filtered data

    Args:
        ripper_folder (str): path to folder containing mzml ripper jsons
        parameters_file (str): full filepath to input parameters json file.
        Defaults to 'C:/Users/group/PolymerMassSpec/Examples/InputParams_Test.json'.
    """
    # create directory for output folder, where filtered ripper_jsons will be
    # saved
    filter_outputdir = os.path.join(ripper_folder, 'filtered_rippers')

    # check whether filter output folder exists; if not, create it
    if not os.path.exists(filter_outputdir):
        os.mkdir(filter_outputdir)

    # retrieve data extrator parameters from parameters file
    extractor_parameters = return_extractor_parameters(parameters_file)
    print(f'extractor_parameters = {extractor_parameters}')

    err = extractor_parameters["error"]
    err_abs = extractor_parameters["err_abs"]

    # get pre_screen filter parameters dict from within extractor parameters
    filter_parameters = extractor_parameters['pre_screen_filters']

    filter_parameters["err"] = err
    filter_parameters["err_abs"] = err_abs

    # get list of mzml ripper jsons in ripper_folder
    ripper_jsons = list_files_by_type(ripper_folder, '.json')

    # iterate through list of ripper_jsons, reading each as a dictionary
    for ripper_json in ripper_jsons:
        ripper_dict = open_json(ripper_json)

        # filter ripper_dict according to pre_screening_filters specified
        # in input parameters .json file
        # ripper_dict = apply_filter_dict_ripper_json(
        #     ripper_dict,
        #     filter_parameters
        # )
        first_spectrum = ripper_dict['ms1']['spectrum_1']
        ripper_dict = {'ms1': first_spectrum}
        # name output file for filtered mass spec data in ripper_dict
        ripper_filename = os.path.basename(ripper_json)
        filtered_json = ripper_filename.replace('.json', '_PRE_FILTERED.json')
        filtered_json = os.path.join(filter_outputdir, filtered_json)

        # finally, write the filtered dict to file
        print('about to write')
        write_to_json(ripper_dict, filtered_json)

def apply_filter_dict_ripper_json(
    ripper_dict,
    filter_dict
):
    """
    Takes a dictionary of mass spec data in mzml_ripper format OR a full
    filepath to an mzml_ripper .json file and returns a filtered mass spec
    data dict (also in mzml_ripper format) based on filters specified in
    filter_dict

    Args:
        ripper_dict (dict or str): mzml_ripper dict or full str filepath to
            mzml_ripper .json file
        filter_dict (dict): dictionary of filter parameters, derived from
            input parameters .json file

    Returns:
        ripper_dict (dict): dictionary of filtered mass spec data, with filters
            from filter_dict applied to data
    """

    # check to see whether ripper_dict is a dictionary or string; if string,
    # it must be a full filepath to an mzml ripper .json file
    if type(ripper_dict) == str:
        ripper_dict = open_json(ripper_dict)

    # extract filter parameters for passing on to later functions
    min_rt = filter_dict['min_rt']
    max_rt = filter_dict['max_rt']

    essential_sigs = filter_dict['essential_signatures']
    signature_ms_level = filter_dict['signature_ms_level']

    # retrieve any ms2 precursors to be applied to filter
    ms2_precursors = filter_dict['ms2_precursors']

    # THIS IS HACKY - CHANGE IT ASAP SO THAT PRECURSORS
    # AND PRECURSOR MS LEVELS ARE READ IN AS DICT
    precursor_ms_level = 2

    # retrieve minimum and maximum MS1, MS2 total intensities
    min_MS1_total_intensity = filter_dict['min_MS1_total_intensity']
    min_MS2_total_intensity = filter_dict['min_MS2_total_intensity']

    # retrieve minimum and maximum MS1, MS2 maximum intensities
    min_MS1_max_intensity = filter_dict['min_MS1_max_intensity']
    min_MS2_max_intensity = filter_dict['min_MS2_max_intensity']

    # retrieve error threshold and whether threshold is in ppm or absolute
    # mass units
    err = filter_dict['err']
    err_abs = filter_dict['err_abs']

    # return filtered ripper_dict with spectra not meeting filter criteria
    # removed
    ripper_dict = apply_pre_screening_filters(
        ripper_dict,
        min_rt,
        max_rt,
        essential_sigs,
        signature_ms_level,
        ms2_precursors,
        precursor_ms_level,
        min_MS1_total_intensity,
        min_MS2_total_intensity,
        min_MS1_max_intensity,
        min_MS2_max_intensity,
        err,
        err_abs
    )

    return ripper_dict
