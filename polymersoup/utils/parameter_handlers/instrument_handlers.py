"""
This file contains constants and functions for setting and retrieving default
parameter values for mass spectrometers and chromatography used to generate
data for PolymerSoup
"""
import os
import json
HERE = os.path.dirname(__file__)
INSTRUMENT_STANDARDS = os.path.join(
    os.path.dirname(HERE), "instrument_standards")
MASS_SPECS = os.path.join(INSTRUMENT_STANDARDS, "mass_spec_models")
CHROMATOGRAPHY = os.path.join(INSTRUMENT_STANDARDS, "chromatography")

#  list of parameters potentially dependenent on mass spec used
MASS_SPEC_INST_PARAMS = [
    "error",
    "err_abs",
    "rt_units",
    "min_ms1_max_intensity",
    "min_ms2_max_intensity",
    "fragmentation",
    "fragment_series",
    "dominant_signature_cap"
]

#  list of parameters potentially dependent on chromatography
CHROMATOGRAPHY_INST_PARAMS = [
    "min_rt",
    "max_rt",
    "resolution_high",
    "resolution_backup",
    "analytes"
]


#  create dict of aliases and default instrument settings for mass specs
MASS_SPEC_CONFIGS = {
    alias.lower().replace(".json", ""): os.path.join(MASS_SPECS, alias)
    for alias in os.listdir(MASS_SPECS)
}


#  create dict of aliases and default instrument settings for chromatography
CHROMATOGRAPHY_CONFIGS = {
    alias.lower().replace(".json", ""): os.path.join(CHROMATOGRAPHY, alias)
    for alias in os.listdir(CHROMATOGRAPHY)
}

def check_instrument_info(params_dict):
    """
    Check default values for parameters that can be defined by mass spec
    and chromatography instruments.

    Args:
        params_dict (dict): input parameters dict supplied from input file

    Returns:
        dict: dict of instrument info with format:
        {
            "mass_spec": {
                param: mass_spec_default
            },
            "chromatography: {
                param: chromatography_default
            }
        }
    """
    instrument_info = {}
    instrument_info["mass_spec"] = retrieve_instrument_alias_info(
        params_dict=params_dict,
        inst_key="instrument",
        configs=MASS_SPEC_CONFIGS)
    instrument_info["chromatography"] = retrieve_instrument_alias_info(
        params_dict=params_dict,
        inst_key="chromatography",
        configs=CHROMATOGRAPHY_CONFIGS)
    return instrument_info

def retrieve_instrument_alias_info(params_dict, inst_key, configs):
    """
    |Takes alias of mass spec or chromatography instrument and retrieves
    default parameters defined in instrument config file

    Args:
        params_dict (dict): dict of parameters and values from input file
        inst_key (str): either "mass_spec" or "chromatography" to specify which
            type of instrument
        configs (dict): dict of available instrument aliases and associated
            defaults

    Raises:
        Exception: raised if instrument aliases does not match pre-configured
            instruments

    Returns:
        dict: dict of parameters and associated default values for instrument
    """
    if inst_key in params_dict:
        if params_dict[inst_key]:
            if params_dict[inst_key] in configs:
                with open(configs[params_dict[inst_key]]) as f:
                    info = json.load(f)
                    return info
            else:
                raise Exception(
                    f'instrument {params_dict[inst_key]} not\n'
                    'configured. Please choose a valid instrument alias or\n'
                    'create a new instrument config file. Valid aliases:\n'
                    f'{configs.keys}')
    return {}
