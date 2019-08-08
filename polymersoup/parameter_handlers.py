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
    extractor_parameters = read_parameters(parameters_file)[
        "extractor_parameters"]

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

def return_polymer_parameters(
    parameters_file
):
    """
    [summary]
    
    Args:
        parameters_file ([type]): [description]
    
    Returns:
        [type]: [description]
    """
    polymer_config_file = open_json(parameters_file)[
        "directories"]["polymer_config_file"]

    polymer_config_dict = open_json(polymer_config_file)

    return polymer_config_dict