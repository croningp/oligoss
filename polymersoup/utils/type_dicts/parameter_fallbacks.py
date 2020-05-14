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
    "instrument_dependent": {
        "mass_spec": None,
        "chromatography": None
    },
    "chromatography": None
}

#  fallbacks for nested silico parameters
SILICO_LIBRARY_FALLBACKS = {
    "isomeric_targets": None,
    "modifications": None,
    "ms1": {
        "min_z": 1,
        "max_z": None,
        "universal_sidechain_modifications": True,
        "universal_terminal_modifications": True
    },
    "ms2": {
        "adducts": None,
        "signatures": None
    }
}

#  fallbacks for nested extractor parameters
EXTRACTOR_FALLBACKS = {
    "pre_screen_filters": None,
    "min_ms1_total_intensity": None,
    "min_ms2_total_intensity": None,
    "min_ms1_max_intensity": 1000,
    "min_ms2_max_intensity": None
}

#  fallbacks for nested postprocess parameters
POSTPROCESS_FALLBACKS = {
    "essential_fragments": None,
    "exclude_fragments": None
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
