"""
This file contains functions for reading parameters from input files and
creating instance of Parameters class for use in PolymerSoup workflows
"""
import json
import copy
from .parameters import Parameters
from ..type_dicts.parameter_fallbacks import (
    ESSENTIAL_CORE_PARAMS, CORE_CLASSES, CORE_PARAM_FALLBACKS)

def generate_parameters(
    params_json,
    param_classes=CORE_CLASSES
):
    """
    Generates instance of Parameters class

    Args:
        params_json (str): full file path of input parameters file
        param_classes (List[str], Optional): list of parameter subclasses.
            Defaults to CORE_CLASSES

    Raises:
        Exception: exception raised for missing essential parameters

    Returns:
        Parameters: instance of Parameters class
    """

    # read params file as dict
    with open(params_json) as f:
        params_dict = json.load(f)

    #  init dict to store formatted parameters dict, ready to be passed into
    #  Parameters
    formatted_params = {}

    #  iterate through parameter inner classes and add to formatted_params
    for param_class in param_classes:
        formatted_params[param_class] = generate_param_subobject(
            params_dict=params_dict,
            param_class=param_class)

    #  make sure core parameters not specific to nested classes are in
    #  formatted_params to be included as properties in final instance of
    #  Parameters class
    for param in ESSENTIAL_CORE_PARAMS:
        if param not in formatted_params:
            formatted_params[param] = params_dict[param]
    for param in CORE_PARAM_FALLBACKS:
        if param not in formatted_params:
            if param in params_dict:
                formatted_params[param] = params_dict[param]
            else:
                formatted_params[param] = CORE_PARAM_FALLBACKS[param]

    for param in CORE_PARAM_FALLBACKS:
        if param not in formatted_params:
            if param in params_dict:
                formatted_params[param] = params_dict[param]
            else:
                formatted_params[param] = CORE_PARAM_FALLBACKS[param]

    #  create instance of Parameters ready to be used in other modules
    #  and return
    return Parameters(
        params_dict=formatted_params,
        params_class="core")

def generate_param_subobject(params_dict, param_class):
    """
    Creates instance of Parameters to be used as an inner class

    Args:
        params_dict (dict): dict of parameters and values
        param_class (str): id of inner parameter class

    Returns:
        Parameters: instance of Parameters class
    """
    # there is a weird bug whereby the params dict is changed after going
    # through parameters, making copies is a temporary fix as it's removing some
    # keys from silico ms1 currently
    params_silico = copy.deepcopy(params_dict)
    params_silico_ms1 = copy.deepcopy(params_dict)
    params_silico_ms2 = copy.deepcopy(params_dict)

    #  get parameter class from inputs. NOTE: silico parameters is the only
    #  parameter class with nested inner classes (ms1 and ms2)
    if param_class != "silico":
        return Parameters(
            params_dict=params_dict,
            params_class=param_class)

    #  generate instance of Parameters class for silico parameters
    silico = Parameters(
        params_dict=params_silico,
        params_class="silico")

    #  set silico.ms1 to instance of Parameters class for silico ms1 parameters
    silico.ms1 = Parameters(
        params_dict=params_silico_ms1,
        params_class="silico_ms1")

    #  set silico.ms2 to instance of Parameters class for silico ms2 parameters
    silico.ms2 = Parameters(
        params_dict=params_silico_ms2,
        params_class="silico_ms2")

    if not silico.ms2.adducts:
        silico.ms2.adducts = silico.ms1.adducts

    #  return silico parameters
    return silico
