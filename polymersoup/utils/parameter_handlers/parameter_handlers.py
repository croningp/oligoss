"""
This file contains functions for reading parameters from input files and
creating instance of Parameters class for use in PolymerSoup workflows
"""

from ..type_dicts.parameter_type_dicts import (
    SILICO_LIBRARY_TYPES,
    EXTRACTOR_TYPES,
    POSTPROCESS_TYPES,
    CORE_PARAM_TYPES
)
from ..type_dicts.parameter_instrument_dependencies import (
    SILICO_LIB_DEPENDENCIES,
    EXTRACTOR_DEPENDENCIES,
    POSTPROCESS_DEPENDENCIES
)
from .instrument_handlers import check_instrument_info
from .parameters import Parameters

#  list of parameters which must ALWAYS be supplied in input file
ESSENTIAL_PARAMS = [
    "mode",
    "screening_method",
    "monomers",
    "silico",
    "extractors",
    "postprocess",
    "polymer_class"
]

def generate_parameters(
    params_dict,
    param_classes
):
    """
    Generates instance of Parameters class

    Args:
        params_dict (dict): dict of parameters and associated values
        param_classes (List[str]): list of parameter subclasses

    Raises:
        Exception: exception raised for missing essential parameters

    Returns:
        Parameters: instance of Parameters class
    """

    #  init dict to store formatted parameters dict, ready to be passed into
    #  Parameters
    formatted_params = {}

    #  check for any essential parameters that are missing from input file
    missing_essentials = [
        param for param in ESSENTIAL_PARAMS
        if param not in params_dict]
    if missing_essentials:
        raise Exception(
            f'{len(missing_essentials)} essential parameters missing'
            f'from input parameters file: {missing_essentials}')

    #  get default values for chosen mass spectrometer and / or chromarography
    instrument_info = check_instrument_info(params_dict)

    #  iterate through parameter inner classes and add to formatted_params
    for param_class in param_classes:
        formatted_params[param_class] = generate_param_subobject(
            params_dict=params_dict,
            param_class=param_class,
            instrument_info=instrument_info)

    #  make sure core parameters not specific to nested classes are in
    #  formatted_params to be included as properties in final instance of
    #  Parameters class
    for param in ESSENTIAL_PARAMS:
        if param not in formatted_params:
            formatted_params[param] = params_dict[param]

    #  create instance of Parameters ready to be used in other modules
    #  and return
    return Parameters(
        params_dict=formatted_params,
        type_dict=CORE_PARAM_TYPES,
        instrument_info=instrument_info,
        instrument_dependencies=None)

def generate_param_subobject(params_dict, param_class, instrument_info):
    """
    Creates instance of Parameters to be used as an inner class

    Args:
        params_dict (dict): dict of parameters and values
        param_class (str): id of inner parameter class
        instrument_info (dict): default values for parameters defined by
            instruments (mass spec, chromatography) used in experiment

    Returns:
        Parameters: instance of Parameters class
    """

    #  get types of each value, and instrument dependencies for filling in
    #  non-specified values for parameters
    type_dict, dependencies = retrieve_param_info(param_class)

    #  get parameter class from inputs. NOTE: silico parameters is the only
    #  parameter class with nested inner classes (ms1 and ms2)
    if param_class != "silico":
        return Parameters(
            params_dict=params_dict[param_class],
            type_dict=type_dict,
            instrument_dependencies=dependencies,
            instrument_info=instrument_info)

    #  add ms1 and ms2 parameters to silico parameters
    ms1_params = Parameters(
        params_dict=params_dict["silico"]["ms1"],
        type_dict=type_dict["ms1"],
        instrument_dependencies=dependencies,
        instrument_info=instrument_info)
    ms2_params = Parameters(
        params_dict=params_dict["silico"]["ms2"],
        type_dict=type_dict["ms2"],
        instrument_dependencies=dependencies,
        instrument_info=instrument_info)

    #  get final silico params dict ready to create instance of Parameters class
    #  for silico parameters
    silico_params = {
        "ms1": ms1_params,
        "ms2": ms2_params
    }

    #  return instance of Parameters for silico parameters
    return Parameters(
        params_dict=silico_params,
        type_dict={"ms1": Parameters, "ms2": Parameters},
        instrument_dependencies=dependencies,
        instrument_info=instrument_info)

def retrieve_param_info(param_class):
    """
    Retrieves type_dict and instrument dependencies for a parameter class
    Args:
        param_class (str): str id of parameter class

    Raises:
        Exception: raises exception if parameter is not pre-defined in
            type_dicts and instrument dependencies

    Returns:
        dict, dict: type_dict, instrument dependencies
    """

    #  remove "_parameters" substring if supplied in input file
    param_class.replace("_parameters", "")

    #  return dependencies from parameter class id
    if param_class.lower() == "silico":
        return SILICO_LIBRARY_TYPES, SILICO_LIB_DEPENDENCIES
    elif param_class.lower() == "extractors":
        return EXTRACTOR_TYPES, EXTRACTOR_DEPENDENCIES
    elif param_class.lower() == "postprocess":
        return POSTPROCESS_TYPES, POSTPROCESS_DEPENDENCIES
    else:
        raise Exception(f'{param_class} is invalid')
