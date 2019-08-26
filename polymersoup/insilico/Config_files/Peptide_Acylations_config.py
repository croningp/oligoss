"""
This is a config file for peptide acylating groups

"""
from ..Constants.GlobalChemicalConstants import *

# MODIFICATIONS_DICT: dictionary of covalent modifications, with corresponding
# three letter codes. This dict should provide all relevant information for
# terminal and side chain-specific covalent modifications that can be attached
# to polymers (e.g. N-terminal peptide acylations)
# MODIFICATIONS_DICT is in format:
# {"Mod": {
#       "mass": mass,
#       "terminus": 0, -1, None or 2,
#       "side chain attachments": [monomer one letter codes],
#       "free_mod_fragments": {
#               "pos": [mass],
#               "neg": [mass]
#       },
#       "mass_diff": {
#               "ms1": MS1_mass_diff,
#               "ms2": MS2_mass_diff
#       },
#       "universal_ms2_shift": bool
#       }
# }
# where:
#  "Mod" = covalent modification three letter code.

#  "mass" = neutral monoisotopic mass of modififying group

#  "terminus" = terminus/termini that can be modified by "Mod" (0 for start
#   terminus, -1 for end terminus, 2 for both or None for neither)

#  "side chain attachments": list of monomer one letter codes for monomers that
#  can have the covalent attachment at their side chain

#  "free_mod_fragments" = dict of MS2 fragment m/z values corresponding to free
#   modification group fragments (e.g. (Ole-OH)+ for free oleic acid in
#   positive mode)

# "mass_diff": mass difference when adding modification to MS1 and MS2 ions.
# final masses will be in format: (M + Mod - mass_diff) where M = mass of
# unmodified sequence, Mod = mass of modification and mass_diff = mass_diff

# "universal_ms2_shift": specifies whether all MS2 fragments for subsequences
# containing the modification MUST have the modification accounted for in their
# mass; this is set to False for most modifications as MS2 fragments can often
# be observed + / - modification mass

"""
MODIFICATIONS_DICT for peptides below; see comment beneath dict for full
description of modifying groups
"""
MODIFICATIONS_DICT = {
    #  "Ole" = covalent modification three letter code for oleic acid
    "Ole": {

    # "mass" = neutral monoisotopic mass of oleic acid
        "mass": 282.25589,

    # "terminus" = terminus acylated by oleic acid; 0 = N-terminus (i.e.
    # starting terminus)
        "termini": [0],

    # "side_chain_attachments" = list of monomer one letter codes for monomers
    # with side chains that can be acylated by oleic acid
        "side_chains_attachments": ["K", "R"],

        # "free_mod_fragments" = m/z values of free oleic acid fragments
        # produced from acylated precursors in both positive and negative mode
        "free_mod_fragments": {
            "pos": [265.2532],
            "neg": [281.2181]
        },

        # "mass_diff" = the mass lost when adding oleic acid to a peptide -
        # water as this is a condensation reaction
        "mass_diff": {
            "ms1": H2O,
            "ms2": H2O
        },
        # "universal_ms2_shift" = specifies whether ALL MS2 fragments with
        # acylated termini MUST have the oleic acid mass shift; set to False
        # as acylated fragments can be found + / - oleic acid
        "universal_ms2_shift": False
    },
    "Pal": {
        "mass": 256.24024,
        "terminus": 0,
        "side chain attachments": ["K", "R"],
        "free_mod_fragments": {
            "pos": [239.2375],
            "neg": [255.2324]
        },
        "mass_diff": {
            "ms1": H2O,
            "ms2": H2O
        },
        "universal_ms2_shift": False
    },
    "Ace": {
        "mass": 60.02114,
        "terminus": 0,
        "side chain attachments": ["K", "R"],
        "free_mod_fragmnts": {
            "pos": [],
            "neg": []
        },
        "mass_diff": {
            "ms1": H2O,
            "ms2": H2O
        },
        "universal_ms2_shift": False
    },
    "Fmc": {
        "mass": 240.07866,
        "terminus": -1,
        "side chain attachments": ["D", "E"],
        "free_mod_fragments": {
            "pos": [],
            "neg": []
        },
        "mass_diff": {
            "ms1": H2O,
            "ms2": H2O
        },
        "universal_ms2_shift": False
    },
    "Tbu": {
        "mass": 74.07317,
        "terminus": None,
        "side chain attachments": [],
        "free_mod_fragments": {
            "pos": [],
            "neg": []
        },
        "mass_diff": {
            "ms1": H2O,
            "ms2": H2O
        },
        "universal_ms2_shift": False
    },
    "Trt": {
        "mass": 260.12012,
        "terminus": None,
        "side chain attachments": ["S", "T"],
        "free_mod_fragments": {
            "pos": [],
            "neg": []
        },
        "mass_diff": {
            "ms1": H2O,
            "ms2": H2O
        },
        "universal_ms2_shift": False
    },
    "Boc": {
        "mass": 118.06301,
        "terminus": 0,
        "side chain attachments": [],
        "free_mod_fragments": {
            "pos": [],
            "neg": []
        },
        "mass_diff": {
            "ms1": H2O,
            "ms2": H2O
        },
        "universal_ms2_shift": False
    }
}
"""
Ole = oleic acid; Pal = palmitic acid; Ace = acetic acid
Fmc = Fmc; Tbu = tert-butyl; Trt = trityl; Boc = boc
"""

# TERMINAL_LINKERS = dictionary of modification groups that can cross-link
# sequences at multiple termini. This dict is used for modifications that have
# more than one reactive point (e.g. organic acids with two carboxylates for
# peptide acylation)

# TERMINAL_LINKERS format= {
#   "Mod": {
#       "mass": mass
#       "termini": {
#           "0": n,
#           "-1": x
#        },
#       "mass_diff": {
#           "ms1": {
#               "0": Mass,
#               "-1": Mass
#           },
#           "ms2": {
#               "0": Mass,
#               "-1": Mass
#           },
#       "universal_ms2_shift": bool
#       }
#   }
# }
# where:
# "Mod" = modification linker three letter code. This should be unique, and
# MUST NOT BE FOUND IN MODIFICATIONS_DICT

# "mass" = neutral monoisotopic mass of linker group

# "termini" = subdict to specify terminal linkages and number of links that can
# be formed with each terminus type. format = {"0": n, "-1": x} where "0" =
# key to specify start terminus of polymer(s); "-1" = key to specify end
# terminus of polymer(s); n = number of links that can be formed with start
# termini (i.e. number of linked termini per Mod); x = number of links that
# can be formed with end termini (i.e. number of linked termini per Mod)

# "mass_diff": mass difference when adding modification to MS1 and MS2 ions.
# final masses will be in format: (M + Mod - mass_diff) where M = mass of
# unmodified sequence, Mod = mass of modification and mass_diff = mass_diff

# "universal_ms2_shift": specifies whether all MS2 fragments for subsequences
# containing the modification MUST have the modification accounted for in their
# mass; this is set to False for most modifications as MS2 fragments can often
# be observed + / - modification mass


TERMINAL_LINKERS = {
    # "Suc" = three letter code for succinic acid, which has two carboxylates
    # and is therefore able to link two sequences by their N-termini
    # (terminus 0)
    "Suc": {

        # "mass" = neutral monoisotopic mass of succinic acid
        "mass": 118.02663,

        # "termini" = subdict to specify which termini can be linked by succinic
        # acid, and how many of each termini it can link; as there are two
        # carboxylates, succinic acid can link two N-termini (terminus 0) and 0
        # C-termini (terminus -1)
        "termini": {
            "0": 2,
            "-1": 0
        },

        # "mass_diff" = mass lost when adding linker modification to termini.
        # This is water for succinic acid, as acylation at the N-terminus is
        # a condensation reaction; -1 terminus mass_diff is set to None
        # as succinic acid should not react with the C-terminus
        "mass_diff": {
            "ms1": {
                "0": H2O,
                "-1": None
            },
            "ms2": {
                "0": H2O,
                "-1": None
            }
        },
        "universal_ms2_shift": False
    }
}
"""
"Suc" = succinic acid
"""
