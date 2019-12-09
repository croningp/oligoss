"""

This file contains standard parameters for extracting and processing data 
from mass spectrometry instruments used commonly in polymer sequencing 
experiments. 

"""

# below is a dictionary of standard data extraction parameters for mass 
# spectrometers used commonly in polymer sequencing experiments 
INSTRUMENT_STANDARDS = {
    
    "Orbitrap": {

        "extractor_parameters": {

            "error": 0.01,
            "err_abs": True,
            "min_ms2_peak_abundance": 90,
            "pre_run_filter": True,
            "min_MS1_total_intensity": None,
            "min_MS2_total_intensity": None,
            "min_MS1_max_intensity": 1E5,
            "min_MS2_max_intensity": None,

            "pre_screen_filters": {
                "min_rt": 0,
                "max_rt": 1E6,
                "essential_signatures": None,
                "signature_types": ["Im"],
                "signature_ms_level": 2,
                "massdiff_bins": False,
                "ms2_precursors": None,
                "min_MS1_total_intensity": None,
                "min_MS2_total_intensity": None,
                "min_MS1_max_intensity": None,
                "min_MS2_max_intensity": None
            }
        }
    },

    "Bruker": {

        "extractor_parameters": {
            "error": 0.01, 
            "err_abs": True,
            "min_ms2_peak_abundance": 90,
            "pre_run_filter": True,
            "min_MS1_total_intensity": None,
            "min_MS2_total_intensity": None,
            "min_MS1_max_intensity": 1E4,
            "min_MS2_max_intensity": 1E2,
        
            "pre_screen_filters": {
                    "min_rt": 0,
                    "max_rt": 1E6,
                    "essential_signatures": None,
                    "signature_types": ["Im"],
                    "signature_ms_level": 2,
                    "massdiff_bins": False,
                    "ms2_precursors": None,
                    "min_MS1_total_intensity": None,
                    "min_MS2_total_intensity": None,
                    "min_MS1_max_intensity": None,
                    "min_MS2_max_intensity": None
            }
        }
    }
}

# below is a list of all data extraction parameters that must be defined 
# for each experiment - either by instrument defaults or custom inputs by the
# user 
DATA_EXTRACTION_PARAMETERS = [
    "error",
    "err_abs",
    "min_ms2_peak_abundance",
    "pre_run_filter",
    "min_MS1_total_intensity", 
    "min_MS2_total_intensity",
    "min_MS1_max_intensity",
    "min_MS2_max_intensity",
    "pre_screen_filters",
    "N_bpc_peaks"
]