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
H = 1.00783
H2O = 18.010565
OH = 17.00274
NH3 = 17.02655

# Peptide Backbone Functional Groups

CO = 27.99492
COOH = 44.99767
NH = 15.0109
NH2 = 16.01872
CH3CN = 41.02655
CH3OH = 32.02622

SOLVENTS = {
    "H2O": H2O, # water
    "ACN": CH3CN, # acetonitrile
    "MeOH": CH3OH # methanol
}

"""
Section 2: Adducts

This section contains information on common cationic and anionic adducts.
CATIONS and ANIONS contain information on positive and negative ions,
respectively. The general format of these dictionaries is as follows:
{f'{symbol}' : (mass, min_z, max_z)} where symbol = atomic symbol or chemical
formula, mass = monoisotopic mass, min_z and max_z = minimum and maximum absolute
charge for a single adduct, respectively

"""
CATIONS = {
    "H": (1.00783, 1, 1),
    "Na": (22.989218, 1, 1),
    "K": (38.963158, 1, 1),
    "NH4": (18.033823, 1, 1),
    "Mn": (54.93805, 2, 7),
    "Ca": (39.96259, 1, 2),
    "Zn": (63.92915, 1, 2),
    "Fe": (53.93961, 1, 3),
    "Co": (58.9332, 2, 3)
}
ANIONS = {
    "Cl": (34.969402, 1, 1),
    "Br": (78.918885, 1, 1),
    "I":  (126.90448, 1, 1),
    "ClO4": (98.94853, 1, 1),
    "-H": (-1.00783, 1, 1)
}
