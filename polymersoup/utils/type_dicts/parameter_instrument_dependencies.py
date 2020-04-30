"""
This file contains parameter class instrument-dependencies
"""

# parameters to be retrieved from instrumentation if not specified for building
# in silico libraries
SILICO_LIB_DEPENDENCIES = {
    "mass_spec": ["fragment_series", "mode", "max_neutral_losses"],
    "chromatography": []
}

# paramters to be retrieved from instrumentation if not specified for building
# screening data
EXTRACTOR_DEPENDENCIES = {
    "mass_spec": [
        "error",
        "min_ms2_peak_abundance",
        "min_ms1_max_intensity",
        "min_ms2_max_intensity",
        {
            "pre_screen_filters": [
                "min_ms1_max_intensity",
                "min_ms2_max_intensity"]}
    ],
    "chromatography": ["min_rt", "max_rt"]
}

# parameters to be retrieved from instrumentation if not specified for post-
# processing screened data
POSTPROCESS_DEPENDENCIES = {
    "mass_spec": [
        "core_linear_series",
        "dominant_signature_cap"],
    "chromatography": [
        "rt_bin",
        "ms2_rt_bin"
    ]
}
