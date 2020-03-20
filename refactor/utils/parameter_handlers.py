import os
import json
from typing import List, Dict, Optional
from global_chemical_constants import FUNCTIONAL_GROUPS
from instrument_handlers import (
    MASS_SPEC_INST_PARAMS,
    CHROMATOGRAPHY_INST_PARAMS,
    MASS_SPEC_CONFIGS,
    CHROMATOGRAPHY_CONFIGS)

silico_type_dict = {
    "mode": str,
    "max_length": int,
    "min_length": int,
    "monomers": List[str],
    "optional": None,

    "MS1": {
        "ms1_adducts": List[str],
        "min_z": int,
        "max_z": int,
        "losses": bool,
        "max_neutral_losses": int,
        "universal_sidechain_modifications": bool,
        "universal_terminal_modifications": bool,
        "terminal_modifications": Optional[Dict[str, str]],
        "side_chain_modifications": Optional[Dict[str, List[str]]],
        "cyclic_sequences": bool,
        "isomeric_targets": Optional[List[str]],
        "optional": ["terminal_modifications", "side_chain_modifications"]
    },

    "MS2": {
        "fragment_series": List[str],
        "ms2_adducts": Optional[List[float]],
        "ms2_losses": bool,
        "add_signatures": bool,
        "signatures": Optional[List[str]],
        "min_z": int,
        "max_z": int,
        "optional": ["ms2_adducts", "signatures"]
    }
}

extractor_type_dict = {
    "error": float,
    "error_units": str,
    "min_ms2_peak_abundance": Optional[float],
    "pre_run_filter": bool,
    "min_ms1_total_intensity": Optional[float],
    "min_ms1_max_intensity": Optional[float],
    "pre_filter": bool,
    "rt_units": str,
    "optional": ["min_ms1_total_intensity", "min_ms1_max_intensity"],

    "pre_screen_filters": {
        "min_rt": Optional[float],
        "max_rt": Optional[float],
        "essential_signatures": Optional[List[float]],
        "signature_ms_level": Optional[int],
        "min_ms1_max_intensity": Optional[float],
        "min_ms1_total_intensity": Optional[float],
        "min_ms2_max_intensity": Optional[float],
        "min_ms2_total_intensity": Optional[float],
        "optional": "all"
    }
}

postprocess_type_dict = {
    "exclude_frags": Optional[List[str]],
    "optional_core_frags": Optional[List[str]],
    "dominant_signature_cap": float,
    "essential_fragments": Optional[List[str]],
    "subsequence_weights": List[float]
}

def load(json_file):
    with open(json_file) as f:
        return json.load(f)

def generate_params_dict(params_file):

    with open(params_file) as f:
        load = json.load(f)
    return load

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
    if (not params_dict["extractor_parameters"]["instruments"]
            or "instruments" not in params_dict["extractor_parameters"]):
        return None, None
    instrument_info = params_dict["extractor_parameters"]["instruments"]
    if "mass_spec" in instrument_info:
        if type(instrument_info["mass_spec"]) == str:
            if (instrument_info["mass_spec"].lower().replace(".json", "")
                    not in MASS_SPEC_CONFIGS):
                if not os.path.isfile(instrument_info["mass_spec"]):
                    raise Exception(f'filepath to mass spectrometer config file\n'
                        f'({instrument_info["mass_spec"]}) is not a valid\n'
                        'filepath. Please enter either an alias for one of\n'
                        'the pre-made mass spectrometer config files or\n'
                        'or provide a valid filepath. Here are a list of pre-\n'
                        'configured aliases for mass spectrometers:\n'
                        f'{MASS_SPEC_CONFIGS.keys()}')
            mass_spec_info = load(instrument_info["mass_spec"])

def generate_extractor_params(params_dict):
    """[summary]

    Args:
        params_dict ([type]): [description]
    """
    extractor_params = params_dict["extractor_parameters"]
    for param in extractor_type_dict:
        if param not in extractor_params:
            if not extractor_params["instrument"]
def retrieve_fg_value(param_value):

    if type(param_value) == str:
        if param_value[0]:
            if param_value[0] == "-":
                return -FUNCTIONAL_GROUPS[param_value]
            return FUNCTIONAL_GROUPS[param_value]

def format_params(param_value, target_type):

    if type(param_value) == target_type:
        return param_value

    if target_type == float or List[float]:
        if type(param_value) == str:
            return retrieve_fg_value(param_value)
        if type(param_value) == List[str]:
            return [
                retrieve_fg_value(param_value=value)
                for value in param_value]
