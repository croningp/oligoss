"""
This file contains constants and functions for setting and retrieving default
parameter values for mass spectrometers and chromatography used to generate
data for PolymerSoup
"""
import os
import json

from .polymer_param_handlers import load_polymer_info
from ..errors import InvalidMSFragmentation

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

def retrieve_instrument_alias_info(params_dict, inst_key, configs=None):
    """
    Takes alias of mass spec or chromatography instrument and retrieves
    default parameters defined in instrument config file

    Args:
        params_dict (dict): dict of parameters and values from input file
        inst_key (str): either "mass_spec" or "chromatography" to specify which
            type of instrument.
        configs (dict, optional): dict of available instrument aliases and
            associated defaults.

    Raises:
        Exception: raised if instrument aliases does not match pre-configured
            instruments

    Returns:
        dict: dict of parameters and associated default values for instrument
    """
    if inst_key in params_dict:
        if not configs:
            if inst_key == "instrument":
                configs = MASS_SPEC_CONFIGS
            elif inst_key == "chromatography":
                configs = CHROMATOGRAPHY_CONFIGS
            else:
                raise Exception(
                    f'{inst_key} not a valid instrument parameter key')

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

def retrieve_mass_spec_info_params_obj(params_obj):

    if not params_obj.instrument:
        return {}

    with open(MASS_SPEC_CONFIGS[params_obj.instrument]) as f:
        info = json.load(f)
        return info

def sanity_check_silico_fragmentation(
    silico_params,
    polymer_class,
    ms_level,
    instrument_info,
    fragmentation_method='default'
):
    """
    Checks that silico fragmentation parameters are compatible with chosen
    polymer_class and mass spec instrumentation - i.e. makes sure all proposed
    fragments could exist in a real mass spec!

    Args:
        silico_params (dict): parameters dict for relevant silico params
        polymer_class (str): polymer class alias
        ms_level (int): ms level to check fragmentation
        instrument_info (dict): dict of mass spec parameters
        fragmentation_method (str, optional): type of fragmentation method that
            is being used. If "default", it is assumed that the fragmentation
            method could be any valid method for instrument.
            Defaults to 'default'.

    Raises:
        InvalidMSFragmentation (Exception): raised if any fragments are
            incompatible with instrument and polymer class in any way
    """

    #  retrieve polymer config file by looking up polymer alias in
    #  polymer_configs folder
    polymer_config = load_polymer_info(polymer_class)

    #  get available fragmentation methods for instrument at chosen ms level
    valid_fragmentations = instrument_info["fragmentation"][f'ms{ms_level}']

    #  confirm that chosen fragmentation method is available in instrument
    if fragmentation_method == "default":
        pass
    elif fragmentation_method not in valid_fragmentations:
        raise Exception(
            f'Fragmentation via {fragmentation_method} is not available\n'
            'on your mass spectrometer. Please choose one of the\n'
            f'following: {valid_fragmentations.keys()}')
    else:
        valid_fragmentations = [fragmentation_method]

    #  check that neutral losses are permitted
    if "neutral" not in valid_fragmentations:
        if silico_params["max_neutral_losses"]:
            raise InvalidMSFragmentation(
                frag_method=valid_fragmentations,
                neutral_invalid=True)
    valid_fragmentations.remove("neutral")

    #  get available fragmentation pathways for chosen polymer class
    polymer_frags = polymer_config["FRAG_SERIES"]

    #  get compatible fragmentation pathways that can induce fragmentations in
    #  chosen polymer class
    permissible_fragment_pathways = polymer_frags["default_linear"]

    #  check that chosen fragment series are available for polymer class
    for series in silico_params["fragment_series"]:
        if series not in polymer_frags:
            raise InvalidMSFragmentation(
                frag_method=valid_fragmentations,
                fragment_series=series)

        #  check that fragmentation method(s) are also compatible with polymer
        for frag_method in valid_fragmentations:
            if frag_method not in permissible_fragment_pathways:
                raise InvalidMSFragmentation(
                    frag_method=frag_method)

            #  check that each fragment series is compatible with the chosen
            #  combination(s) of polymer and fragmentation method
            if series not in permissible_fragment_pathways[frag_method]:
                raise InvalidMSFragmentation(
                    frag_method=frag_method,
                    fragment_series=series)
