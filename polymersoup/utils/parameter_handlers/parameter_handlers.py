"""
This file contains functions for reading parameters from input files and
creating instance of Parameters class for use in PolymerSoup workflows
"""
from .parameters import Parameters
from ..type_dicts.parameter_fallbacks import ESSENTIAL_CORE_PARAMS, CORE_CLASSES

def generate_parameters(
    params_dict,
    param_classes=CORE_CLASSES
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

    #  get parameter class from inputs. NOTE: silico parameters is the only
    #  parameter class with nested inner classes (ms1 and ms2)
    if param_class != "silico":
        return Parameters(
            params_dict=params_dict,
            params_class=param_class)

    #  return instance of Parameters for silico parameters
    return Parameters(
        params_dict=params_dict,
        params_class="silico")
