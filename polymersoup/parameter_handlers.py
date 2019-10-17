"""
This file contains functions for extracting info on run parameters for
in silico fragmentation, data extraction and post-processing.

"""
import time

# import filehandler - the file containing functions that load parameters
# from files
from .filehandler import *

# import data extraction and postprocessing standards for mass spec instruments
# commonly used in experiments
from .StandardInstrumentParameters.InstrumentStandards import *

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
    """
    Reads essential, monomer-specific signature ions (denoted by monomer one-
    letter codes) from input parameters and returns list of m/z values for 
    essential signature ions. 
    
    Args:
        parameters_file (str): full file path to input paremeters file
    
    Returns:
        list: list of m/z values for signature ions
    """
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
    parameters = read_parameters(parameters_file)
    extractor_parameters = parameters["extractor_parameters"]

    # check for default instrument and create dict of default data extraction
    # parameters for that instrument
    if "instrument" in parameters.keys():
        instrument = parameters['instrument']
        instrument_params = INSTRUMENT_STANDARDS[instrument][
            'extractor_parameters']
    else:
        instrument = None

    extractor_params = {}

    # check that every data extraction parameter is defined in parameters file
    # if not, attempt to load it from default instrument parameters
    # this this is not possible, prompt the user to provide the parameter values
    for parameter in DATA_EXTRACTION_PARAMETERS: 
        try: 
            extractor_params[parameter] = extractor_parameters[parameter]
        except KeyError: 
            
            if instrument:
                print(f'parameter={parameter}')
                print(f'instrument_params={instrument_params.keys()}')
                extractor_params[parameter] = instrument_params[parameter]
            
            else:
                print(f'you are missing a value for {parameter} for the {instrument} mass spectrometer')
                parameter_value = input('Please enter its value here:')
                extractor_params[parameter] = parameter_value

    # update extractor params with polymer-specifc signature ions used
    # in extracting data and sequence screening
    extractor_params[
        'pre_screen_filters']['essential_signatures'] = load_essential_signature_ions(
            parameters_file
        )          

    return extractor_params

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

    # add date and time signature to parameters dict 
    parameters_dict = add_timestamp()

    # combine all above into a final parameters_dict that can be read by
    # the executable
    parameters_dict.update(
        {
            "silico_parameters": silico_dict,
            "extractor_parameters": extractors_dict,
            "postprocessing_parameters": postprocess_dict,
            "directories": directories_dict,
            "screening_method": read_parameters["screening_method"],
            "data_extraction": read_parameters["data_extraction"],
            "postprocess": read_parameters["postprocess"]
        }
    )

    # return all parameters in dict
    return parameters_dict

def add_timestamp(
    day=True,
    month=True,
    year=True
):
    """
    Returns date and time dict for saving to run_parameters log files
    
    Args:
        day (bool, optional): specify whether to log day. Defaults to True.
        month (bool, optional): specify whether to log month. Defaults to True.
        year (bool, optional): specify whether to log year. Defaults to True.
    
    Returns:
        dict: dict of date and time. s
    """
    timestamp_dict = {}

    date_string = ''

    if day: 
        date_string += f'{time.strftime("%d")}/'
    
    if month:
        date_string += f'{time.strftime("%m")}/'
    
    if year:
        date_string += f'{time.strftime("%y")}'
    
    timestamp_dict["date"] = date_string

    hour = round(float(time.strftime("%H")), 2)
    minute = round(float(time.strftime("%M")), 0)
    second = round(float(time.strftime("%S")), 0)
    

    timestamp_dict["time"] = f'{hour}h:{minute}m:{second}s'

    return timestamp_dict
