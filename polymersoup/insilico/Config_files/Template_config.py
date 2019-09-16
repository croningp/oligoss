"""
This file is a template for a polymer constants file.
Once filled in, the information contained herein should contain all constants
needed to define MS and MSMS fragmentation properties of chosen polymer class
to build MS and MSMS sequence libraries for screening mass spectrometry data

"""

from .. import GlobalChemicalConstants


"""
SECTION 1: GENERAL MS1 AND CHEMISTRY RULES
This section should contain the basic information for constructing in silico
MS1 data for polymer class.

MONOMERS = dictionary of monomer one letter codes, along with associated neutral
            monoisotopic masses and reactivity classes (functional groups).
            This dictionary should specify the neutral monoisotopic mass for
            each monomer, its reactive functional groups, and the number of each
            reactive functional group. This information is essential to ensure
            that all theoretical sequences generated are chemically feasible.
            Format =
                    {
                        'X': [
                            mass, [
                                [rxn_class1, n],
                                [rxn_class2, y]
                            ]
                    }
            where X = monomer one letter code; mass = neutral monoisotopic
            mass; rxn_class1 and rxn_class2 = reactivity classes (functional
            groups - e.g. amine, aldehyde); n = number of rxn_class1 functional
            groups; y = number of rxn_class2 functional groups.


MASS_DIFF = the mass difference when adding an additional monomer on to a polymer
            chain (e.g. for condensation polymer = WATER)

ELONGATION_UNIT = the number of additional monomer units typically added when
                elongating a polymer. This will typically be 1 for most polymers

REACTIVITY_CLASSES = dictionary of reactivity classes with associated compatible
                classes and monomers. Format = {'classA' : [['classX', 'classY'],
                                                ['A', 'B']}
                where classA = reaction class, classX and classY = classes
                that are cross-reactive with classA, and A and B are monomers
                within reaction class classA

SYMMETRY = bool to define whether polymer is identical at both ends. Set to false
            for polymers with different termini (e.g. N- and C- termini for peptides),
            true for polymers with identical functional groups at both termini
            in linear chains

CHAIN_TERMINATORS = list of monomers that terminate chain elongation

LOSS_PRODUCTS = dictinary of monomer one letter codes and associated side chain
            neutral loss products, i.e. masses that can be lost from the monomer
            side chain.

IONIZABLE_SIDECHAINS = dictionary of monomers that can be ionised with extra adducts
    at the SIDE CHAIN, with associated adducts, minimum and maximum absolute charge
    states in both positive and negative mode. Format = {"X": {"pos": (adduct, a, b),
    "neg": (adduct, a, b)}} where "X" = monomer one letter code, adduct = adduct string
    (must be found in either CATIONS OR ANIONS in GlobalChemicalConstantss), a = min 
    side chain charge, b = max side chain charge for ionized form.

INTRINSICALLY_CHARGED_MONOMERS = dictionary of monomers that have an intrinsic charge
    (i.e. charged without addition of adducts), with associated lists of permissible adducts
"""

MONOMERS = {"X": [0, [["rxn_class1", 1], ["rxn_class2", 1]]]}

"""
INSERT DETAILS REGARDING MONOMERS HERE - FULL NAMES, CAS NUMBERS (IF APPLICABLE),
ADDITIONAL RELEVANT DETAILS

"""

MASS_DIFF: float = 0
ELONGATION_UNIT = 1
REACTIVITY_CLASSES = {
    "rxn_class": [["rxn_class1", "rxn_class2"], ["X", "Y"]]
}
"""
INSERT ADDITIONAL DETAILS REGARDING REACTIVITY CLASSES HERE - FUNCTIONAL GROUPS
AND ANY OTHER RELEVANT DETAILS

"""

SYMMETRY = False
"""
INSERT ADDITIONAL DETAILS REGARDING POLYMER TERMINI HERE - FUNCTIONAL GROUPS,
ETC
"""

CHAIN_TERMINATORS = ["C", "D"]

"""
INSERT ADDITIONAL DETAILS REGARDING CHAIN TERMINATORS HERE

"""

LOSS_PRODUCTS = {
    "X": []
}

IONIZABLE_SIDECHAINS = {
    "X": {"pos": (),
          "neg": ()}
}

INTRINSICALLY_CHARGED_MONOMERS = {}
"""
IONIZABLE_SIDECHAINS: Side chains with free amines (K, R, H) can gain a proton
and become cationic; side chains with free carboxylates (E, D) can lose a proton
and become anionic
"""
"""
INSERT ADDITIONAL DETAILS REGARDING SIDE CHAIN LOSS PRODUCTS HERE

"""

"""
SECTION 2: MS2 FRAGMENTATION

This section should contain all the information required to construct a basic
MS2 fragment series for linear polymers. Each fragment type is defined as a key
in the FRAG_SERIES dict. Typically, fragment keys are one letter codes used to
denote fragments - e.g. the standard 'b' and 'y' fragment series for peptides.
When building fragment series, the fragment generator will use the following
convention:
                f'{frag}{n}'
frag = fragment one letter code, n = number of monomers in fragment (e.g. b1,
        b2, b3, y1, y2, y3 etc...)

FRAG_SERIES = dictionary of fragment series one letter codes and associated
properties. Each fragment series key = fragment one letter code, value =
subdictionary containing all relevant fragment properties. Below is a list of
fragment properties defined in the FRAG_SERIES dict. This information has been
filled in for the well studied peptide b- and y- fragment series as an example.

Fragment properties:

'terminus' {0 or -1} - specifies which terminus the fragment series starts from.
                    0 = start of sequence, -1 = end of sequence.

'mass_diff' {float} - the difference in mass between the fragment and corresponding
                    subsequence. This may vary depending on whether cations or
                    anions are being fragmented, therefore separate mass_diffs
                    are specified for positive mode ('pos') and negative mode
                    ('neg') mass spec. Example: the 'y3' fragment of peptide
                    sequence 'AGVS' = the mass of (GVS+H) in positive mode, and
                    (GVS-H) in negative mode.

'fragmentation_unit' - the minimum number of monomer units typically added and
                    / or removed at a time when building a fragment series.
                    Default = ELONGATION_UNIT (see Section 1)

'start' - the starting position of the fragment series within the fragmenting
            sequence, relative to the first monomer in the sequence. If this is
            0, fragment series will be generated along the full length of the
            sequence; otherwise it will be generated from 0+n position, where
            n = 'start'

'end' - the ending position of the fragment series within the fragmenting
            sequence, relative to the last monomer in the sequence. If this is
            0, fragment series will be generated along the full length of the
            sequence; otherwise it will be generated up until -(n+1) position,
            where n = 'end'

'intrinsic_adducts' - OPTIONAL. this subdictionary contains any intrinsic
            adducts generated as default for a given fragment series in
            positive and negative mode. Intrinsic adducts are adducts which
            are added to the fragment as default (for example: +H for peptide y
            fragments in positive mode). This information is important when
            adding non-standard adducts (e.g. Na+, K+, Cl- etc..) to MS2
            fragments, as intrinsic adduct masses must be removed when adding
            other adduct masses

'permissible_adducts' - OPTIONAL. this subdictionary contains lists of adducts
            that may be found associated with a particular MS2 fragment series.
            This is to be used only in cases where there are restrictions on
            associated fragments, as otherwise it is assumed that any MS2
            adducts included in the in silico run are compatible with all
            fragment series. Example: peptide b fragments are inherently
            charged as acylium ions, and as they are extremely unlikely to be
            multiply charged in the absence of ionizable or intrinsically
            charged side chains, there are no permissible adducts for this
            series.

MS2_SIGNATURE_IONS = MS2 fragments which can be used as markers for monomers
                    and / or small subsequences.
"""


FRAG_SERIES = {
    "b": {
        "terminus": 0,
        "mass_diff": {
            "pos": -OH,
            "neg": -H
        },
        "fragmentation_unit": {
            "pos": ELONGATION_UNIT,
            "neg": ELONGATION_UNIT
        },
    "y": {
        "terminus": -1,
        "mass_diff": {
            "pos": H,
            "neg": -H
        },
        "fragmentation_unit": {
            "pos": ELONGATION_UNIT,
            "neg": ELONGATION_UNIT
        }
    }
}
}

MS2_SIGNATURE_IONS = {}

"""

"""
