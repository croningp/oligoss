.. _OLIG-Files:

##########
OLIG Files
##########

OLIG files define the full scope of known ionization and fragmentation pathways for an oligomer class.

Layout
======

A proper OLIG configuration file should contain all information required to generate full possible sequence libraries for target oligomer class, as well as MS1 precursors and MS2 product ions.
It should be divided into three sections:

#. **General Oligomer Properties**:

    This defines properties required for generating MS1 precursor libraries.

#. **MS2 Fragmentation Properties**:

    This defines properties required for generating MS2 product ions.

#. **Modifications**:

    This defines any covalent modifications and their appropriate targets.

.. _General-Properties:
General Oligomer Properties
==========================

#. **MONOMERS**:
    * Description: defines monomer one-letter codes, their associated monoisotopic neutral masses and reactive functional groups.
    * Type: `Dict[str, list]`
    * Options: N/A
    * Example: {"A": [89.04768, [["amine", 1], ["carboxyl", 1]], "alanine"]} defines the amino acid monomer alanine, with a neutral monoisotopic mass of 89.04768, 1 reactive amine and 1 reactive carboxylic acid.

#. **MASS_DIFF**:
    * Description: defines mass lost (or gained) upon addition of a monomer to an elongating chain.
    * Type: `float` or `str`
    * Options: any valid float corresponding to neutral monoisotopic mass difference in standard mass units (u) *or* a functional group string (e.g. "H2O") corresponding to the mass difference.

#. **ELONGATION**:
    * Description: defines standard number of monomer additions per elongation event (in monomer units).
    * Type: `int`
    * Options: any valid integer >= 1.

#. **REACTIVITY_CLASSES**:
    * Description: defines cross-reactivity of monomer functional groups.
    * Type: `Dict[str, List[list]]`
    * Options: keys must correspond to functional groups found in monomer library.
    * Example: The following example defines cross-reactivity of the "amine" functional group for an oligomer class, which in this case can react with either "carboxyl" or "hydroxyl" groups. The "amine" functional group is found in monomers "A", "B" and "C": {"amine": ["carboxyl", "hydroxyl"], ["A", "B", "C"]}.

#. **SYMMETRY**:
    * Description: determines whether termini are functionally equivalent (i.e. whether forward sequence == reverse sequence).
    * Type: `bool`
    * Options: true or false.
    * Example: **SYMMETRY** = true for peptides as they have distinct C- and N-termini.
#. **LOSS_PRODUCTS**:
    * Description: specifies any neutral loss fragmentation events that can occur for specific monomer sidechains.
    * Type: `Dict[str, list]`
    * Options: keys must correspond to valid monomer one-letter codes, values must be lists of either monoisotopic neutral mass losses (`float`) or functional group strings corresponding to these neutral losses.
    * Example: {"N": ["NH3", "H2O"]} for monomer "N" with possible neutral loss fragmentations corresponding to loss of either an ammonia ("NH3") or water ("H2O") mass.

  .. note::
    It is assumed that these fragmentation events can either be a product of in-source CID at MS1 or standard MS2 CID.


#. **IONIZABLE_SIDECHAINS**:
    * Description: specifies non-backbone ionization sites that occur at specific monomer sidechains.
    * Type: `Dict[str, dict]`
    * Options: keys must correspond to valid monomer one-letter codes. Values must be dictionaries with keys "pos" and "neg" defining possible ionization events for positive and negative mode, respectively.
    * Example: In the following example, monomer "K" can be ionized via proton addition in positive mode and the monomer "D" can be ionized via proton abstraction in negative mode: {"K": {"pos": ["H", 1, 1], "neg": null}, "D": {"pos": null,"neg": ["-H", 1, 1]}}

#. **INTRINSICALLY_CHARGED_MONOMERS**:
    * Description: this defines monomers which have an intrinsic, non-exchangeable charge.
    * Type: `Dict[str, int]`
    * Options: keys must be valid monomer one-letter codes, with values equivalent to intrinsic charge state due to non-exchangeable ions.
    * Example: for a monomer "Z" with intrinsic charge of -2: {"Z": -2}

#. **SIDE_CHAIN_CROSSLINKS**:
    * Description: this defines any monomer-monomer crosslinks that can occur, and their effects on MS1 ionization and MS2 fragmentation.
    * Type: `Dict[str, dict]`
    * Options: keys must correspond to valid monomer one-letter codes. Key-Value pairs in the subdict are as follows:

      * *monomers*:
           * Description: defines other monomers that can form sidechain crosslinks with target monomer.
           * Type: `List[str]`
           * Options: list of valid monomer one-letter codes.

      * *crosslink_massdiff*:
           * Description: defines mass lost or gained upon crosslinking.
           * Type: `float` or `str`
           * Options: either a float corresponding to neutral monoisotopic mass diff or a string representing functional group mass diff.

      * *permissible_crosslink_charges*:
           * Description: defines permissible charge states for crosslinked moiety at the sidechain(s) of crosslinked monomers.
           * Type: `List[int]`
           * Options: list of any valid integer corresponding to permissible charge states.

      * *disrupt_ms2*:
           * Description: specifies whether crosslinking event disrupts standard linear fragmentation along backbone.
           * Type: `bool`
           * Options: true or false.

    * Example: The following example defines cross-linking events for the monomer "K" which, in its non-crosslinked state can be ionized at its sidechain (see **IONIZABLE_SIDECHAINS**). It can crosslink with monomers "E" and "D" via sidechain links. However, this type of crosslinking event does not disrupt standard linear MS2 fragmentation pathways: {"K": {"monomers": ["E", "D"], "crosslink_massdiff": "H2O", "permissible_crosslink_charges": [0], "disrupt_ms2": false}}

.. MS2-Properties:

MS2 Fragmentation Properties
============================

MS2 fragmentation properties are required for defining possible MS2 fragmentation pathways for an oligomer class.
This includes both linear fragment series and signature ion fragments.

#. **FRAG_SERIES**:
    * Description: dictionary that defines all linear fragmentation pathways. Linear fragmentation pathways are defined as any fragment series indexed stepwise along the oligomer backbone.
    * Type: `Dict[str, dict]`
    * Options: keys must correspond to fragment series one-letter codes. Properties of individual fragment series are defined in subparameters (see **Defining FRAG_SERIES**, below).
    * Example: see **Defining  FRAG_SERIES** section.
#. **MS2_SIGNATURE_IONS**:
    * Description: this defines any monomer-specific signature ions that may occur. 
    * Type: `Dict[str, list]`
    * Options: Keys correspond to signature ion str code, values lists of monomer one-letter codes and corresponding free signature m/z values.
    * Example: {"Im": ["F", 120.0813], ...} defines "Im" signature fragment for monomer "F" with m/z 120.0813.
    
    .. note::
        The same fragmentation events that produce linear fragment series can also produce signature ions.
        However, Polymersoup considers these as separate events due to the diversity of possible signature ions.

#. **MODIFICATIONS**:
    * Description: defines any covalent modifications and possible modification sites.
    * Type: `Dict[str, dict]`
    * Options: Keys must be strings corresponding to modification three-letter codes. Values are subdicts defining modification properties (see **Defining MODIFICATIONS**, below).

.. _Frag-Series:

Defining FRAG_SERIES
--------------------

The **FRAG_SERIES** dict is used to define properties relevant to linear fragment series (i.e. fragment series that are indexed stepwise along the oligomer backbone).

#. **default_linear**:
    * Description: specifies default linear fragmentation pathways to be included in silico libraries depending on the mass spec fragmentation method used to acquire data.
    * Type: `Dict[str, List[str]]`
    * Options: keys must correspond to valid fragmentation methods defined in **Instrument_Configs**. Values are list of linear fragment series codes for linear fragment series that are produced via the specified fragmentation method.
    * Example: In the case of "HCD" fragmentation producing "a", "b" and "y" MS2 fragments: {"HCD": ["b", "y", "a"]}.

    .. note::
        * There is redundancy with **Instrument Configs**.
        * Linear fragment series can also be specified for individual oligomer classes in instrument config files.
        * These can also be overwritten directly in input parameters file.

#. **default_core**:
    * Description: specifies default core linear fragment series (i.e. linear series used in confidence assignments) depending on the mass spec fragmentation method used to acquire data.
    * Type: `Dict[str, List[str]]`
    * Options: keys must correspond to valid fragmentation methods defined in **Instrument Configs**. Values are lists of linear fragment series codes for linear fragment series that are produced via the specified fragmentation method *and* are required for assigning confidence scores.
    * Example: In the case of "HCD" fragmentation producing core fragment series "b" and "y": {"HCD": ["b", "y"]}.
    
    .. note::
        * There is redundancy with **Instrument Configs**.
        * Linear fragment series can also be specified for individual oligomer classes in instrument config files.
        * These can also be overwritten directly in input parameters file.

#. **terminus**:
    * Description: specifies "home" terminus from which linear fragment series is indexed.
    * Type: `int`
    * Options: either 0 or -1 for fragment series indexed from terminus 0 and -1, respectively.

    .. note::
        For oligomer classes with **symmetry** == False, terminus is irrelevant.

#. **mass_diff**:
    * Description: specifies *neutral* mass difference between a fragment and its corresponding intact neutral sequence slice.
    * Type: `str` or `float`
    * Options: valid float corresponding to mass difference in mass units (u) or functional group string representing a neutral monoisotopic mass corresponding to mass (e.g. "OH", "H2O").

#. **fragmentation_unit**:
    * Description: specifies increment of fragment indices when producing linear fragment series.
    * Type: `Dict[str, [int or str]`
    * Options: must be a key for "pos", "neg" defining fragmentation unit in positive and negative mode, respectively. Values must either be ints >= 0 or strings representing int value (most commonly "ELONGATION_UNIT" if **fragmentation_unit** == **ELONGATION_UNIT**).
    * Example: **fragmentation_unit** will equal 1 or **ELONGATION_UNIT** for the majority of oligomer classes. A possible exceptions to this would be for alternating copolymers with alternating backbone links:

#. **start**:
    * Description: start position of fragment series relative to **terminus**.
    * Type: `int`
    * Options any valid integer >= 0.
    * Example: 0 for a fragment series that begins immediately at **terminus**, 1, 2, 3 for fragment series that begins 1, 2 or 3 indices away from **terminus**.

#. **end**:
    * Description: end position of fragment series relative to other terminus (i.e. terminus 0 and -1 for **terminus** == 1 and **terminus** == 0, respectively).
    * Type: `int`
    * Options: any valid integer >= 0. 
    * Example: 0 for a fragment series that terminates at final index on backbone, 1, 2, 3 for fragment series that terminates 1, 2 or 3 indexes away from final index on backbone.

#. **intrinsic_charge**:
    * Description: defines any *non-exchangeable* ions associated with fragments of a particular series.
    * Type: `Dict[str, int]`
    * Options: keys must be "pos" and "neg" for positive and negative mode, respectively. Values must be integers representing intrinsic charge value.
    * Example: {"pos": 1, "neg": null} for a fragment series with intrinsic charge of 1 in positive mode but no intrinsic charge in negative mode.

    .. note::
        Do not confuse this with **intrinsic_adduct**. By definition MS2 fragment series are charged by default. However, this can be a result of either *non-exchangeable* or *exchangeable* ions. **intrinsic_charge** defines charge state due to *non-exchangeable* ions.

    .. note::
        **NOT_REQUIRED**. This property does not need to be defined if **intrinsic_adduct** is defined. However, at least one of these properties must be defined to account for fragment charge.

#. **intrinsic_adduct**:
    * Description: defines any *exchangeable* ions associated with fragments of a particular series.
    * Type: `Dict[str, str]`
    * Options: keys must be "pos" and "neg" for positive and negative mode, respectively. Values must correspond to strings representing adducts stored in **Global Chemical Constants**.
    * Example: {"pos": "H", "neg": "-H"} for a fragment series that is intrinsically protonated in positve mode but deprotonated in negative mode.
    
    .. note::
        This property should only be used to define *exchangeable ions* (i.e. ions that can be swapped for extrinsic ions in sample matrix). Do not confuse with *non-exchangeable ions_, which are defined in  **intrinsic_charge**.

#. **exceptions**:
    * Descriptions: for oligomer classes with mixed backbones (i.e. more than one backbone bond type that can be fragmented at MS2), fragmentation properties may differ depending on what type of bond is being fragmented at a particular index.
    * Type: `Dict[str, dict]`
    * Options: keys must be "pos" and "neg" to define exceptions to standard fragmentation rules in positive and negative mode, respectively. Values define exceptions to any combination of previously described MS2 fragmentation properties for linear fragment series.
    * Format: {mode (str): {func_group: {prop: {"positions": List[int], "start": int, "end": int, "exception_value: Value}}}}.

      * *mode* == either "pos" or "neg" for positive or negative mode.
      * *func_group* == functional group that causes exception to standard fragmentation pathway.
      * *prop* == the property for which the exception may apply.
      * *positions*: defines list of indexes in subsequence at which exception applies. Some fragmentation exceptions only apply when the non-standard backbone link is in a particular position in the fragment subsequence.
      * *start*: defines start position at which exception applies, relative to home terminus.
      * *end*: specifies number of indices away from end terminus at which exception no longer applies
      * *exception_value*: the substituted value to use for the property if exception applies.
    * Example: The following example is for a fragment series with exception to **mass_diff** in cases where a bond between a "hydroxyA"-containing monomer is being fragmented. The exception applies when the "hydroxyA"-containing monomer occurs at the final index of the subsequence. The exception applies from the very first index of the fragment series but ends one index away from the end terminus: {"pos": {"hydroxyA": {"mass_diff": {"positions": [-1], "start": 0, "end": 1, "exception_value": 26.98709}}}}
