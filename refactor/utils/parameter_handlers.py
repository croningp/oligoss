import os
import json
from typing import List
from global_chemical_constants import FUNCTIONAL_GROUPS
from instrument_handlers import (
    MASS_SPEC_CONFIGS,
    CHROMATOGRAPHY_CONFIGS)
from type_dicts.parameter_type_dicts import (
    silico_type_dict,
    extractor_type_dict,
    postprocess_type_dict)

def generate_params_dict(params_file):

    with open(params_file) as f:
        load = json.load(f)
    silico_params = generate_silico_params(load)
    return silico_params

def generate_silico_params(params_dict):
    """
    Generates and formats parameters relevant to building in silico libraries
    for a screening experiment
    Args:
        params_dict (dict): dict read from user input parameters file

    Raises:
        Exception:
        Exception:

    Returns:
        silico_params (dict): dict of all parameters needed to build in silico
            libraries for sequencing
    """
    silico_params = params_dict["silico_parameters"]

    for param, param_value in silico_params.items():
        if param not in ["MS1", "MS2"]:
            if not param_value and param not in silico_params["optional"]:
                raise Exception(f'value for {param} is missing from silico\n'
                                f'parameters. This parameter should be of\n'
                                f'type: {silico_type_dict[param]}')
            if param_value:
                param_value = format_params(
                    param_value=param_value,
                    target_type=silico_type_dict[param])
        else:
            for sub_param, sub_value in param_value.items():
                if not sub_value and sub_param not in silico_params[
                        sub_param]["optional"]:
                    raise Exception(f'value for {sub_param} is missing from\n'
                                    f'{param} silico parameters. This\n'
                                    f'parameter should be of type:\n'
                                    f'{silico_type_dict[param][sub_param]}')
                if sub_value:
                    sub_value = format_params(
                        param_value=sub_value,
                        target_type=silico_type_dict[sub_param][sub_value])
    return silico_params

def load_instrumentation(params_dict):
    """[summary]

    Args:
        params_dict ([type]): [description]

    Raises:
        Exception: [description]
        Exception: [description]

    Returns:
        [type]: [description]
    """
    if (not params_dict["extractor_parameters"]["instruments"]
            or "instruments" not in params_dict["extractor_parameters"]):
        return None, None
    instrument_info = params_dict["extractor_parameters"]["instruments"]
    if "mass_spec" in instrument_info:
        if type(instrument_info["mass_spec"]) == str:
            if (instrument_info["mass_spec"].lower().replace(".json", "")
                    not in MASS_SPEC_CONFIGS):
                if not os.path.isfile(instrument_info["mass_spec"]):
                    raise Exception(f'filepath to mass spectrometer config\n'
                                    f'({instrument_info["mass_spec"]}) is not\n'
                                    f'a valid filepath. Please enter an alias\n'
                                    f'for one of the pre-made mass spect\n'
                                    f'config files or provide a valid\n'
                                    f'filepath. Here are a list of aliases\n'
                                    f'for pre-configured instruments:\n'
                                    f'{MASS_SPEC_CONFIGS.keys()}')
            # mass_spec_info = load(
                # MASS_SPEC_CONFIGS[instrument_info["mass_spec"]])
        if "chromatography" in instrument_info:
            if (instrument_info["chromatography"].lower().replace(".json", "")
                    not in CHROMATOGRAPHY_CONFIGS):
                if not os.path.isfile(instrument_info["mass_spec"]):
                    raise Exception(f'file path to chromatography config file\n'
                                    f'is invalid. Please enter a valid\n'
                                    f'alias or a filepath to a valid\n'
                                    f'chromatography config file. Here are a\n'
                                    f'list of aliases for pre-configured\n'
                                    f'chromatography separations:\n'
                                    f'{CHROMATOGRAPHY_CONFIGS.keys()}')

def retrieve_fg_value(param_value):
    """[summary]

    Args:
        param_value ([type]): [description]

    Returns:
        [type]: [description]
    """
    if type(param_value) == str:
        if param_value[0]:
            if param_value[0] == "-":
                return -FUNCTIONAL_GROUPS[param_value]
            return FUNCTIONAL_GROUPS[param_value]

def format_params(param_value, target_type):
    """[summary]

    Args:
        param_value ([type]): [description]
        target_type ([type]): [description]

    Returns:
        [type]: [description]
    """
    if type(param_value) == target_type:
        return param_value

    if target_type == float or List[float]:
        if type(param_value) == str:
            return retrieve_fg_value(param_value)
        if type(param_value) == List[str]:
            return [
                retrieve_fg_value(param_value=value)
                for value in param_value]
