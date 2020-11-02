"""
This file contains global chemical constants - i.e. masses and other associated,
unchanging properties for various elements and molecules encountered frequently
in the ALife's mass spectrometry related experiments

"""

"""
Section 1: Common Non-Adduct Masses

This section contains monoisotopic masses of important molecules and ions that
are not adducts.
SOLVENTS = dictionary of comment neutral solvent molecules and associated
masses

"""
# Quick constants that will be used often as neutral losses and mass_diffs
H = 1.00783
H2O = 18.010565
OH = 17.00274
NH3 = 17.02655
CO = 27.99492
COOH = 44.99767
NH = 15.0109
NH2 = 16.01872
CH3CN = 41.02655
CH3OH = 32.02622
C2H5OH = 46.04187
CO2 = 43.98984

# Solvent Dictionary
SOLVENTS = [
    H2O,  # water
    CH3CN,  # acetonitrile
    CH3OH,  # methanol
    C2H5OH  # ethanol
]

# Functional Group Dictionary
FUNCTIONAL_GROUPS = {
    "H2O": H2O,  # water
    "h2o": H2O,  # water
    "water": H2O,  # water
    "NH3": NH3,  # ammonia
    "H": H,  # proton
    "OH": OH,  # hydroxyl group
    "COOH": COOH,  # carboxylic acid group
    "NH2": NH2,  # primary amino group
    "NH": NH,  # secondary amino group
    "CO": CO,  # carbonyl / CO
    "CH3CN": CH3CN,  # acetonitrile
    "MeCN": CH3CN,  # acetonitrile
    "mecn": CH3CN,  # acetonitrile
    "ACN": CH3CN,  # acetonitrile
    "acetonitrile": CH3CN,  # acetonitrile
    "MeOH": CH3OH,  # methanol
    "meoh": CH3OH,  # methanol
    "methanol": CH3OH,  # methanol
    "CH3OH": CH3OH,  # methanol
    "EtOH": C2H5OH,  # ethanol
    "etoh": C2H5OH,   # ethanol
    "ethanol": C2H5OH,  # ethanol
    "CO2": CO2,  # carbon dioxide
    "co2": CO2   # carbon dioxide
}

"""
Section 2: Adducts

This section contains information on common cationic and anionic adducts.
CATIONS and ANIONS contain information on positive and negative ions,
respectively. The general format of these dictionaries is as follows:
{f'{symbol}' : (mass, min_z, max_z)} where symbol = atomic symbol or chemical
formula, mass = monoisotopic mass, min_z and max_z = minimum and maximum
absolute charge for a single adduct, respectively

"""
CATIONS = {
    "H": (1.00783, 1, 1),
    "Na": (22.989218, 1, 1),
    "K": (38.963158, 1, 1),
    "NH4": (18.033823, 1, 1),
    "Cu": (62.9296, 1, 3),
    "Mn": (54.93805, 2, 7),
    "Ca": (39.96259, 1, 2),
    "Zn": (63.92915, 2, 3),
    "Fe": (55.9349375, 2, 3, 4, 6),
    "Co": (58.9332, 2, 4),
    "Ni": (57.93535, 2, 3)
}

ANIONS = {
    "Cl": (34.969402, 1, 1),
    "Br": (78.918885, 1, 1),
    "I": (126.90448, 1, 1),
    "ClO4": (98.94853, 1, 1),
    "-H": (-1.00783, 1, 1),
}
