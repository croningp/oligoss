from typing import List, Dict, Optional

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