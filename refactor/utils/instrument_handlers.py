import sys, os
HERE = os.path.dirname(__file__)
mass_specs = os.path.join(HERE, ".", "instrument_standards", "mass_specs")
chromatography = os.path.join(
    HERE,
    ".",
    "instrument_standards",
    "chromatography")

MASS_SPEC_INST_PARAMS = [
    "error",
    "err_abs",
    "rt_units",
    "min_MS1_max_intensity",
    "min_MS2_max_intensity"
]

CHROMATOGRAPHY_INST_PARAMS = [
    "min_rt",
    "max_rt",
    "resolution_high",
    "resolution_backup"
]

MASS_SPEC_CONFIGS = {
    alias: os.path.join(mass_specs, alias)
    for alias in os.listdir(mass_specs)
}

CHROMATOGRAPHY_CONFIGS = {
    alias: os.path.join(chromatography, alias)
    for alias in os.listdir(chromatography)
}
