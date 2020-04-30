"""
This file contains type_dicts for parameters. Values for parameters provided
in input files must either match specified types in the appropriate type_dict
or be convertible to this type
"""

from typing import List, Dict
from ..parameter_handlers.parameters import Parameters

#  types for parameters in core Parameter class
CORE_PARAM_TYPES = {
    "mode": str,
    "monomers": List[str],
    "silico": Parameters,
    "extractors": Parameters,
    "postprocess": Parameters,
    "screening_method": str,
    "polymer_class": str
}


#  types for nested silico parameters
SILICO_LIBRARY_TYPES = {
    "max_length": int,
    "min_length": int,
    "optional": None,

    "ms1": {
        "adducts": List[str],
        "min_z": int,
        "max_z": int,
        "losses": bool,
        "max_neutral_losses": int,
        "universal_sidechain_modifications": bool,
        "universal_terminal_modifications": bool,
        "terminal_modifications": Dict[str, str],
        "side_chain_modifications": Dict[str, List[str]],
        "cyclic_sequences": bool,
        "isomeric_targets": List[str],
        "optional": ["terminal_modifications", "side_chain_modifications"]
    },

    "ms2": {
        "fragment_series": List[str],
        "adducts": List[float],
        "max_neutral_losses": bool,
        "add_signatures": bool,
        "signatures": List[str],
        "min_z": int,
        "max_z": int,
        "optional": ["ms2_adducts", "signatures"]
    }
}


#  types for nested extractor parameters
EXTRACTOR_TYPES = {
    "error": float,
    "error_units": str,
    "min_ms2_peak_abundance": float,
    "pre_run_filter": bool,
    "min_ms1_total_intensity": float,
    "min_ms1_max_intensity": float,
    "min_ms2_total_intensity": float,
    "min_ms2_max_intensity": float,
    "pre_filter": bool,
    "rt_units": str,
    "optional": [
        "min_ms1_total_intensity",
        "min_ms1_max_intensity",
        "min_ms2_total_intensity",
        "min_ms2_max_intensity",
        "pre_screen_filters"],

    "pre_screen_filters": Dict[str, float]
}


#  types for nested postprocess
POSTPROCESS_TYPES = {
    "exclude_frags": List[str],
    "optional_core_frags": List[str],
    "dominant_signature_cap": float,
    "essential_fragments": List[str],
    "subsequence_weight": List[float],
    "core_linear_series": List[str],
    "rt_bin": float,
    "ms2_rt_bin": float,
    "backup_rt_bin": float,
    "ms2_rt_flexible": bool,
    "min_viable_confidence": float
}
