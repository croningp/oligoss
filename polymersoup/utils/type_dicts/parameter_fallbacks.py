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
    }
}

#  fallbacks for nested silico parameters
SILICO_LIBRARY_FALLBACKS = {
    "isomeric_targets": None,
    "ms1": {
        "min_z": 1,
        "max_z": None,
        "universal_sidechain_modifications": True,
        "universal_terminal_modifications": True,
        "sidechain_modifications": None,
        "terminal_modifications": None
    },
    "ms2": {
        "min_z": 1,
        "max_z": None,
        "adducts": None
    }
}

#  fallbacks for nested extractor parameters
EXTRACTOR_FALLBACKS = None

#  fallbacks for nested postprocess parameters
POSTPROCESS_FALLBACKS = None

# list of parameters which must ALWAYS be supplied in input file
ESSENTIAL_CORE_PARAMS = [
    "mode",
    "monomers",
    "silico",
    "extractors",
    "screening_method",
    "postprocess",
    "polymer_class"
]

#  ids for parameter subobjects
CORE_CLASSES = [
    "silico",
    "extractors",
    "postprocess"
]
