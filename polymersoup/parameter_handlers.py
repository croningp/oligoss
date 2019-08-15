"""
This file contains functions for extracting info on run parameters for
in silico fragmentation, data extraction and post-processing.

"""

# import filehandler - the file containing functions that load parameters
# from files
from .filehandler import *

def __init__():
    pass

def return_silico_parameters(
    parameters_file
):
    """
    Takes an input parameters file (.json) and returns dictionary of parameters
    relevant to in silico functions

    Args:
        parameters_file (str): full file path to input parameters json fil

    Returns:
        silico_parameters (dict): dictionary of parameters and associated
                    info for in silico functions
    """

    silico_parameters = read_parameters(parameters_file)["silico_parameters"]

    return silico_parameters

def load_essential_signature_ions(
    parameters_file
):

    parameters_dict = read_parameters(parameters_file)

    extractor_params = parameters_dict['extractor_parameters']

    filter_dict = extractor_params['pre_screen_filters']

    essential_signatures = filter_dict['essential_signatures']

    essential_signature_ions = []

    if essential_signatures:
        for essential_signature in essential_signatures:
            for sig_type in filter_dict['signature_types']:
                essential_signature_ions.extend(
                    [
                        float(ion)
                        for ion in MS2_SIGNATURE_IONS[
                        sig_type][essential_signature]
                    ]
                )

    return essential_signature_ions


def return_extractor_parameters(
    parameters_file
):
    """
    Takes an input parameters file (.json) and returns dictionary of parameters
    relevant to data extraction

    Args:
        parameters_file (str): full file path to input parameters json

    Returns:
        extractor_parameters (dict): dictionary of parameters and associated
                    info for data extraction
    """
    # retrieve extractor parameters from input parameters .json file
    extractor_parameters = read_parameters(parameters_file)[
        "extractor_parameters"]

    # update extractor parameters with polymer-specifc signature ions used
    # in extracting data and sequence screening
    extractor_parameters[
        'pre_screen_filters']['essential_signatures'] = load_essential_signature_ions(
            parameters_file
        )
    return extractor_parameters

def return_postprocess_parameters(
    parameters_file
):
    """
    Takes an input parameters file (.json) and returns dictionary of parameters
    relevant to sequence confirmation and post-processing

    Args:
        parameters_file (str): full file path to input parameters json

    Returns:
        postprocess_parameters (dict): dictionary of parameters and associated
                    info for sequence confirmation and postprocessing
    """
    postprocess_parameters = read_parameters(parameters_file)[
        "postprocessing_parameters"]

    return postprocess_parameters

def return_directories(
    parameters_file
):
    """
    Takes an input parameters file (.json) and returns dictionary of key
    file directories essential for executing de novo sequencing experiments

    Args:
        parameters_file (str): full file path to input parameters .json file

    Returns:
        directories (dict): dictionary of important directories required to
            run executable
    """
    directories = read_parameters(parameters_file)["directories"]

    return directories

def generate_parameters_dict(
    parameters_file
):
    """
    Takes an input parameters .json file and returns a dictionary of all
    relevant parameters contained within the file for successfully running
    a full de novo sequencing experiment in executable

    Args:
        parameters_file (str): full file path to input parameters .json file
    Returns:
        parameters_dict (dict): dictionary of all parameters relevant for
            executing de novo sequencing experiments
    """
    # load parameters relevant for in silico operations
    silico_dict = return_silico_parameters(parameters_file)

    # load parameters relevant for data extraction (and filtering)
    extractors_dict = return_extractor_parameters(parameters_file)

    # load parameters relevant for sequence postprocessing
    postprocess_dict = return_postprocess_parameters(parameters_file)

    # load important directories essential for executing experiments
    directories_dict = return_directories(parameters_file)

    read_parameters = open_json(parameters_file)

    # combine all above into a final parameters_dict that can be read by
    # the executable
    parameters_dict = {
        "silico_parameters": silico_dict,
        "extractor_parameters": extractors_dict,
        "postprocessing_parameters": postprocess_dict,
        "directories": directories_dict,
        "screening_method": read_parameters["screening_method"]
    }

    # return all parameters in dict
    return parameters_dict
