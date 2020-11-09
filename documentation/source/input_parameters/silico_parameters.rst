.. _Silico-Parameters:

Silico Parameters
=================

Silico parameters define the properties required to construst theoretical MS1 precursors and MS2 product ions.
There are general silico parameters, as well as parameters specific to MS1 and MS2 ions, defined in silico.ms1 and silico.ms2 subparameters, respectively.

#. **min_length**:
    * Description: specifies minimum length (in monomer units) of oligomers in product mixtures.
    * Type: `int`
    * Options: Any integer value.
    * Default: no Default available
      
   .. note::
      * Length of target oligomers affects size of the sequence space to screen.
      * Limits on sequence space for screening will depend on oligomer class, types of monomers, length of oligomers and computing resources avaialable.

#. **max_length**:
    * Description: specifies maximum length (in monomer units) of oligomers in product mixtures.
    * Type: `int`
    * Options:Any integer value.
    * Default: no Default available

   .. note::
      *  length of target oligomers affects size of the sequence space to screen.
      * Limits on sequence space for screening will depend on oligomer class,    * Types of monomers, length of oligomers and computing resources avaialable.

#. **isomeric_targets**:
    * Description: specifies isomeric sequence pool for targetting. If specified, only sequences isomeric to one or more isomeric targets will be included in the screen.
    * Type: `List[str]`
    * Options: list of any sequence strings that would be generated from input monomers and length distribution.
    * Default: None (**not required**).

#. **modifications**:
    * Description: specifies targets for any covalent modifications.
    * Type: `Dict[str or int, List[str]]`
    * Options: keys specify modification targets (terminal and sidechain).
       * Terminal Keys: -1 or "-1" / 0 or "0" for terminus -1 and 0, respectively.
       * Sidechain Keys: keys must be one-letter codes for monomers with compatible sidechains (see **Polymer Configs**).
       * Modification Values: list of modification strings. Strings must correspond to valid modification alias (see **Polymer Configs**).
    * Example: {"K": ["Ole", "Pal"], 0: ["Ole"]}
    * Default: None (**not required**).

#. **ms1**:
    * Description: specifies parameters for MS1 silico ions.
    * Type `dict`
    * Options: many sub-parameters (see :ref:`**MS1 Silico Parameters**<MS1-Silico>`).

#. **ms2**:
    * Description: specifies parameters for MS2 silico ions.
    * Type `dict`
    * Options: many sub-parameters (see **MS2 Silico Parameters**).

.. _MS1-Silico:

MS1 Silico Parameters
---------------------

MS1 silico parameters define properties required for generating theoretical MS1 precursor ions.

#. **min_z**:
    * Description: specifies minimum **absolute** charge of MS1 precursor ions.
    * Type: `int`
    * Options: Any integer value.
    * Default: Default value specified by instrument. If no instrument is specified, Default value == **1**.

#. **max_z**:
    * Description: specifies maximum **absolute** charge of MS1 precursor ions.
    * Type: `int`
    * Options: Any integer value >= **min_z**.
    * Default: None (**not required**). Instrument defaults can be specified in instrument config file.

#. **universal_sidechain_modifications**:
    * Description: specifies whether all sidechains targeted for modification will be modified by one or more modifying agent.
    * Type: `bool`
    * Options: true or false.
    * Default: Default can be specified in polymer config file (see **Polymer Configs**). If not default is specified in config, default value == **True**.

#. **universal_terminal_modifications**:
    * Description: specifies whether all termini targeted for modification will be modified by one or more modifying agent.
    * Type: `bool`
    * Options: true or false.
    * Default: Default can be specified in polymer config file (see **Polymer Configs**). If not default is specified in config, default value == **True**.

#. **max_neutral_losses**:
    * Description: specifies cap on maximum number of neutral loss fragmentation events at MS1.
    * Type: `int`
    * Options: any integer value.
    * Default: Default can be specified in instrument config file (see **Instrument Configs**). If default is not specified in config, default value == **None**.

#. **adducts**:
    * Description: specifies extrinsic ions present in analyte matrix that will affect MS1 ionization.
    * Type: `List[str]`
    * Options: list of ion strings.
    * Default: varies with instrument and polymer class. If no instrumentor polymer-specific defaults, this must be specified.

   .. note::
      * All ions must be defined in **Global Chemical Constants**.

MS2 Silico Parameters
---------------------

MS2 silico parameters define properties required for generating theoretical MS2 product ions.

#. **fragment series**:
    * Description: linear fragment types to be included in MS2 silico libraries.
    * Type: `List[str]`
    * Options: list of any valid MS2 fragment one letter codes in polymer config file (see **Polymer Configs**).
    * Defaults: Defaults vary depending on specific instrument and polymer class. If no default is specified in instrument config, this **MUST** be supplied in input parameters.

#. **max_neutral_losses**:
    * Description: specifies cap on maximum number of neutral loss fragmentation events at MS2.
    * Type: `int`
    * Options: any integer value.
    * Default: default can be specified in instrument config file (see **Instrument Configs**). If default is not specified in config, default value == **None**.

#. **signatures**:
    * Description: specifies signature ion types to be included in MS2 silico libraries.
    * Type: `List[str]`
    * Options: list of any valid signature ion types in polymer config file.
    * Defaults: Defaults vary depending on instrument and polymer class. If no default is specified in config files, default value == **None**.

#. **min_z**:
    * Description: specifies minimum **absolute** charge of MS2 product ions.
    * Type: `int`
    * Options: Any integer value.
    * Default: Default value specified by instrument. If no instrument is specified, no fallback default value is available so this must be supplied.
