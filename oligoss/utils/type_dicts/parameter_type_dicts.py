"""
This file contains type_dicts for parameters. Values for parameters provided
in input files must either match specified types in the appropriate type_dict
or be convertible to this type
"""

from typing import List, Dict

#  keys in type_dicts which are not used to directly set attributes in instance
#  of Parameters class
NON_ATTR_KEYS = [
    "instrument_only",
    "instrument_polymer",
    "instrument_dependent",
    "optional",
    "core"
]

#  types for parameters in core Parameter class
CORE_PARAM_TYPES = {
    "mode": str,
    "monomers": List[str],
    "silico": "parameters",
    "extractors": "parameters",
    "postprocess": "parameters",
    "screening_method": str,
    "polymer_class": str,
    "instrument": str,
    "chromatography": str,
    "free_cores": int,
    "max_cores": int,
    "data_folder": str,
    "output_folder": str,
    "hard_memory_limit": float,
    "relative_memory_limit": float,
    "instrument_dependent": {
        "instrument_only": {
            "mass_spec": None,
            "chromatography": None
        },
        "instrument_polymer": {
            "mass_spec": None,
            "chromatography": None
        }
    },
    "optional": [
        "instrument",
        "chromatography",
        "max_cores",
        "free_cores",
        "data_folder",
        "output_folder",
        "hard_memory_limit",
        "relative_memory_limit"]
}

#  types for nested silico parameters
SILICO_LIBRARY_TYPES = {
    "max_length": int,
    "min_length": int,
    "isomeric_targets": List[str],
    "optional": ["isomeric_targets"],
    "modifications": Dict[str, List[str]],
    "ms1": {
        "adducts": List[str],
        "min_z": int,
        "max_z": int,
        "max_neutral_losses": int,
        "universal_sidechain_modifications": bool,
        "universal_terminal_modifications": bool,
        "optional": [
            "max_neutral_losses",
            "min_z",
            "max_z",
            "universal_sidechain_modifications",
            "universal_terminal_modifications"],

        "instrument_dependent": {
            "instrument_only": {
                "mass_spec": None,
                "chromatography": None
            },
            "instrument_polymer": {
                "mass_spec": ["max_neutral_losses", "min_z", "max_z"],
                "chromatography": None
            }
        }
    },

    "ms2": {
        "fragment_series": List[str],
        "adducts": List[str],
        "max_neutral_losses": bool,
        "signatures": List[str],
        "min_z": int,
        "max_z": int,
        "optional": [
            "signatures",
            "adducts",
            "fragment_series",
            "min_z",
            "max_z",
            "signatures"
        ],
        "instrument_dependent": {
            "instrument_only": {
                "mass_spec": None,
                "chromatography": None
            },
            "instrument_polymer": {
                "mass_spec": [
                    "fragment_series",
                    "max_neutral_losses",
                    "signatures",
                    "min_z",
                    "max_z"
                ],
                "chromatography": None
            }
        }
    }
}

#  types for nested extractor parameters
EXTRACTOR_TYPES = {
    "error": float,
    "error_units": str,
    "min_rt": float,
    "max_rt": float,
    "rt_units": str,
    "min_ms2_peak_abundance": float,
    "min_ms1_total_intensity": float,
    "min_ms1_max_intensity": float,
    "min_ms2_total_intensity": float,
    "min_ms2_max_intensity": float,
    "optional": [
        "min_ms1_total_intensity",
        "min_ms1_max_intensity",
        "min_ms2_total_intensity",
        "min_ms2_max_intensity",
        "pre_screen_filters"],
    "instrument_dependent": {
        "instrument_only": {
            "mass_spec": [
                "error",
                "min_ms1_total_intensity",
                "min_ms1_max_intensity",
                "min_ms2_total_intensity",
                "min_ms2_max_intensity",
                "rt_units"
            ],
            "chromatography": None
        },
        "instrument_polymer": {
            "mass_spec": ["min_ms2_peak_abundance"],
            "chromatography": None
        }
    },
    "pre_screen_filters": Dict[str, float]
}

# types for nested postprocess parameters
POSTPROCESS_TYPES = {
    "exclude_fragments": List[str],
    "optional_core_fragments": List[str],
    "dominant_signature_cap": float,
    "essential_fragments": List[str],
    "subsequence_weight": float,
    "core_linear_series": List[str],
    "rt_bin": float,
    "ms2_rt_bin": float,
    "spectral_assignment_plots": bool,
    "min_plot_confidence": int,
    "molecular_assembly": dict,
    "optional": [
        "exclude_fragments",
        "optional_core_fragments",
        "essential_fragments",
        "subsequence_weight",
        "rt_bin",
        "ms2_rt_bin",
        "dominant_signature_cap",
        "core_linear_series",
        "spectral_assignment_plots",
        "min_plot_confidence",
        "molecular_assembly"
    ],
    "instrument_dependent": {
        "instrument_only": {
            "mass_spec": None,
            "chromatography": ["rt_bin", "ms2_rt_bin"]
        },
        "instrument_polymer": {
            "mass_spec": [
                "exclude_frags",
                "optional_core_fragments",
                "dominant_signature_cap",
                "core_linear_series"
            ],
            "chromatography": None
        }
    }
}
