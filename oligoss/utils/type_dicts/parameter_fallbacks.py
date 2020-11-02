"""
This file contains default values for parameters when not supplied by
instrument configurations or directly in the input file. Parameters not included
in fallbacks MUST be set by input file or instrument defaults.
"""

#  fallbacks for parameters in core Parameter class
CORE_PARAM_FALLBACKS = {
    "silico": None,
    "extractors": None,
    "postprocess": None,
    "free_cores": 0,
    "max_cores": 12,
    "data_folder": None,
    "output_folder": None,
    "instrument_dependent": {
        "mass_spec": None,
        "chromatography": None
    },
    "chromatography": None,
    "hard_memory_limit": None,
    "relative_memory_limit": 0.90
}

#  fallbacks for nested silico parameters
SILICO_LIBRARY_FALLBACKS = {
    "isomeric_targets": None,
    "modifications": None,
    "ms1": {
        "min_z": 1,
        "max_z": None,
        "universal_sidechain_modifications": True,
        "universal_terminal_modifications": True,
        "max_neutral_losses": None
    },
    "ms2": {
        "adducts": None,
        "signatures": None,
        "max_neutral_losses": None
    }
}

#  fallbacks for nested extractor parameters
EXTRACTOR_FALLBACKS = {
    "pre_screen_filters": None,
    "min_ms1_total_intensity": None,
    "min_ms2_total_intensity": None,
    "min_ms1_max_intensity": 1E3,
    "min_ms2_max_intensity": None
}

#  fallbacks for nested postprocess parameters
POSTPROCESS_FALLBACKS = {
    "essential_fragments": None,
    "exclude_fragments": None,
    "spectral_assignment_plots": False,
    "min_plot_confidence": 70,
    "molecular_assembly": {
        "min_confidence": 70,
        "consensus": True,
        "min_peak_identity": 0.7,
        "combine_precursors": False,
        "ppm_window": 5
    }
}

# list of parameters which must ALWAYS be supplied in input file
ESSENTIAL_CORE_PARAMS = [
    "mode",
    "monomers",
    "silico",
    "extractors",
    "screening_method",
    "postprocess",
    "polymer_class",
    "instrument"
]

#  ids for parameter subobjects
CORE_CLASSES = [
    "silico",
    "extractors",
    "postprocess"
]
