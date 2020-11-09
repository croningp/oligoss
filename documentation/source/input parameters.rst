.. _Input-Parameters:

################
Input Parameters 
################

.. image:: ../img/typical_polymersoup_execution.png
    :width: 600
    :align: center

Required Parameters
===================

Input parameter files must be in JSON format.
The following parameters are required upon execution of every sequencing experiment:

#. **mode**:
    *  Description: specifies the charge of ions.
    *  Type: `str`
    *  Options: either "pos" or "neg" for positive and negative ion mode, respectively.


#. **monomers**:
    * Description: monomers used in experiment.
    * Type: `List[str]`
    * Options: any combination of monomer one letter codes from polymer-specific configuration files.


#. **screening method**:
    * Description: specifies which sequencing workflow to use.
    * Type: `str`
    * Options: currently only option is "exhaustive".


#. **polymer_class**:
    * Description: specifies the backbone of oligomers being fragmented.
    * Type: `str`
    * Options: any valid alias associated with a polymer config file (see **Polymer Configs**).


#. **silico**:
    * Description: parameters that define *in silico* properties.
    * Type: `dict`
    * Options: many sub-parameters must be defined (see **Silico Parameters**).


#. **extractors**:
    * Description: parameters that define properties essential to matching and filtering in MS/MS data.
    * Type: `dict`
    * Options: many sub-parameters must be defined (see **Extractor Parameters**).


#. **postprocess**:
    * Description: parameters that define properties relevant to assigning confidence scores to sequences based on extracted data.
    * Type: `dict`
    * Options: many sub-parameters must be defined (see **Postprocess Parameters**).

   .. note::
      *  In addition to these parameters, optional parameters for mass spec model and chromatography can be defined.
      *  If instruments are not defined, additional parameters must be defined in **silico**, **extractors** and **postprocess** (see *Instrument Configs*).


#. **instrument**:
    * Description: specifies the mass spec model used to acquire MS/MS data.
    * Type: `str`
    * Options: any valid alias for a pre-configured instrument file.


#. **chromatography**:
    * Description: specifies the chromatography used to separate products prior to detection via MS.
    * Type: `str`
    * Options: any valid alias for a pre-configured chromatography file.


Silico Parameters
=================

Silico parameters define the properties required to construst theoretical MS1 precursors and MS2 product ions.
There are general silico parameters, as well as parameters specific to MS1 and MS2 ions, defined in silico.ms1 and silico.ms2 subparamters, respectively.

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
    *  Default: no Default available


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
    * Options: many sub-parameters (see **MS1 Silico Parameters**).

#. **ms2**:
    * Description: specifies parameters for MS2 silico ions.
    * Type `dict`
    * Options: many sub-parameters (see **MS2 Silico Parameters**).

MS1 Silico Parameters
=====================

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
=====================

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

Extractor parameters define properties required for screening observed MS and MS/MS data for theoretical silico ions, and also for filtering MS data prior to screening.

#. **error**:
    * Description: specifies error threshold for matching theoretical ion *m/z* values to observed ions in MS/MS data.
    * Type: `float`
    * Options: Any float value between 0 and 1.
    * Default: Default value specified by instrument config.

   .. note::
      * All ions must be defined in **Global Chemical Constants**.
      * Default: no fallback. Default is available must be specified in input parameters or instrument config.

#. **error_units**:
    * Description: specifies units of error tolerance value.
    * Type: `str`
    * Options: "ppm" for relative error tolerance in parts per million, or "abs" for absolute error tolerance in mass units (u).
    * Default: Default value specified by instrument config.

   .. note::
      *  no fallback default is available must be specified in input parameters or instrument config.

#. **min_rt**:
    * Description: specifies minimum retention time (in minutes) for data to be used for screening.
    * Type: `float`
    * Options: any valid float between 0 and length of acquisition run time.
    * Defaults: Default value specified by chromatography config. If no chromatography is specified, default value == None (**not required**).

#. **min_rt**:
    * Description: specifies maximum retention time (in minutes) for data to be used for screening.
    * Type: `float`
    * Options: any valid float between 0 and length of acquisition run time.
    * Defaults: Default value specified by chromatography config. If no chromatography is specified, default value == None (**not required**).

#. **min_ms1_total_intensity**:
    * Description: specifies minimum total intensity (taken as the sum of intensities in an EIC) for MS1 precursor ions to be confirmed.
    * Type: `float`
    * Options: Any float.
    * Defaults:    * Default value is specified by instrument config. If no instrument is specified, default value == None (**not required**).

#. **min_ms2_total_intensity**:
    * Description: specifies minimum total intensity (taken as the sum of intensities in an EIC) for MS2 product ions to be confirmed.
    * Type: `float`
    * Options: Any float >= 0.
    * Defaults: Default value is specified by instrument config. If no instrument is specified, default value == None (**not required**).

#. **min_ms1_max_intensity**:
    * Description: specifies minimum peak in intensity (taken as the most intensity signal in an EIC) for MS1 precursor ions to be confirmed.
    * Type: float
    * Options: Any float >= 0.
    * Default: Default value is specified by instrument config. If no instrument is specified, default value == None (**not required**).

#. **min_ms2_max_intensity**:
    * Description: specifies minimum peak in intensity (taken as the most intensity signal in an EIC) for MS2 product ions to be confirmed.
    * Type: float
    * Options: Any float >= 0.
    * Default: Default value is specified by instrument config. If no instrument is specified, default value == None (**not required**).

#. **rt_units**:
    * Description: retention time units in mzML files.
    * Type: `str`
    * Options: either "min" or "sec" for minutes and seconds, respectively.
    * Defaults: This value depends on the mass spec manufacturer, and should generally only be used in instrument config files.
   
   .. warning::
      Please do not define this value in input parameters if possible.

   .. note::
      If you are unsure about the retention time units of your mass spec vendor, convert your mzML files to mzmlripper format directly and check the recorded retention times of the first and final spectra.

#. **pre_screen_filters**:
    * Description: specifies retention time and intensity thresholds for pre-filtering ripper data before screening.
    * Type: `Dict[str, float]`
    * Options: several subparameters to define (see **Pre-Screen Filters**, below).
    * Default: None (**not required**).

Pre-Screen Filters
========================

Pre-screen filters are used to remove irrelevant data from raw spectra before screening.

#. **min_ms1_max_intensity**:
    * Description: specifies minimum peak in intensity in raw MS1 spectra for spectra to be included in later screening.
    * Type: `float`
    * Options: any valid float >= 0.
    * Defaults: None (**not required**).

#. **min_ms2_max_intensity**:
    * Description: specifies minimum peak in intensity in raw MS2 spectra for spectra to be included in later screening.
    * Type: `float`
    * Options: any valid float >= 0.
    * Default: None (**not required**).

#. **min_ms1_total_intensity**:
    * Description: specifies minimum total intensity in raw MS1 spectra for spectra to be included in later screening.
    * Type: `float`
    * Options: any valid float >= 0.
    * Default: None (**not required**).

#. **min_ms2_total_intensity**:
    * Description: specifies minimum total intensity in raw MS2 spectra for spectra to be included in later screening.
    * Type: `float`
    * Options: any valid float >= 0.
    * Default: None (**not required**).

#. **min_rt**:
    * Description: specifies minimum retention time for spectra to be included in later screening.
    * Type: `float`
    * Options: any valid float >= 0.
    * Default: None (**not required**).

#. **max_rt**:
    * Description: specifies maximum retention time for spectra to be included in later screening.
    * Type: `float`
    * Options: any valid float >= 0.
    * Default: None (**not required**).

#. **min_ms2_peak_abundance**:
    * Description: specifies the minimum relative abundance for most intense sequence in an MS2 spectrum. If the relative intensity of the most intense match is less than **min_ms2_peak_abundance**, any fragments found in this spectrum will be discarded.
    * Type: `float`
    * Options: any float in range 0-100.

Postprocess Parameters
======================

Postprocess parameters are used to define properties relevant for assigning sequence confidence based on observed data and how it compares to silico data, and also for any postprocessing such as molecular assembly calculations and spectrum plots.

#. **core_linear_series**:
    * Description: core linear fragment types that are to be used to assign confidence for confirmed sequences from MS2 data.
    * Type: `List[str]`
    * Options: list of any valid fragment types in polymer config file.

   .. note::
      * All fragment types in **core_linear_series** must also be present in **silico.ms2.fragment_series**.
      * However, not all fragments in **silico.ms2.fragment_series** have to be in **core_linear_series**.
      * Fragment types specified *in silico* but not in **core_linear_series** will still be screened and recorded, but not used to assign confidence.

#. **exclude_fragments**:
    * Description: specifies fragments that are to be excluded from confidence calculations, even if they have been confirmed.
    * Type: `List[str]`
    * Options: List of any valid fragment ids in **core_linear_series** (see below).
    * Defaults: None (**not required**).

#. **optional_core_fragments**:
    * Description: core fragments that are to be used in confidence calculations _only_ if confirmed. 
    * If absent from extracted data, these fragments will be ignored in confidence assignments and final confidence score will not be lower as a result of their absence. However, if present, they will be used to calculate confidence.
    * Type: `List[str]`
    * Options: list of any fragment ids from core fragments.
    * Default: Default depends on instrument and polymer class. If no default is specified in instrument config, default value == None (**not required**).

#. **essential_fragments**:
    * Description: fragments that are essential for assigning a sequence's confidence score > 0.
    * Type: `List[str]`
    * Options: list of any valid fragment ids from core linear fragments.
    * Default: Default depends on instrument and polymer config. If not specified in configs, default value == None (**not required**).

#. **dominant_signature_cap**:
    * Description: specifies maximum confidence score for sequences missing one or more expected "dominant" signature ion at MS2.
    * Type: `float`
    * Options: Any float in range 0-100.
    * Default: Default value varies by instrument and polymer class. If not specified in configs, default value == 0 (**not required**).

#. **subsequence_weight**:
    * Description: weighting value for mean continuous fragment coverage in final confidence score for confirmed sequences.
    * Type: `float`
    * Options: any float in range 0-1. 0 = confidence based entirely on % confirmed fragments, 1 = confidence based entirely on mean continuous fragment coverage.
    * Default: No Default. This must be specified in input parameters file.

#. **rt_bin**:
    * Description: minimum resolution (in minutes) between peaks in MS1 EICs.
    * Type: `float`
    * Options: Any valid float between 0 and length of acquisition time.

   .. note::
      * This is no longer used. It will be useful when revisiting quantification of individual sequences from EICs.

#. **ms2_rt_bin**:
    * Description: minimum resolution (in minutes) between peaks in MS2 EICs.
    * Type: `float`
    * Options: Any valid float between 0 and length of acquisition time.

   .. note::
      * This is no longer used. It will be useful when revisiting quantification of individual sequences from EICs.

#. **spectral_assignment_plots**:
    * Description: specifies whether to plot annotated MS2 spectra for confirmed sequences.#
    * Type: `bool`
    * Options: true or false. If true, plots will be saved as PNG images in output directory.
    * Default: false (**not required**).

#. **min_plot_confidence**:
    * Description: specifies minimum confidence score for a sequence's spectral assignments to be plotted.
    * Type: `float`
    * Options: any valid float in range 0-100.
    * Default: 70.

#. **molecular_assembly**:
    * Description: specifies parameters for calculating molecular assembly values for confirmed sequences.
    * Type: `dict`
    * Options: several subparameters (see **Molecular Assembly Parameters**).
    * Default: several Defaults for subparameters.

Molecular Assembly Parameters
=============================

Defines properties relevant for calculating Molecular Assembly for confirmed sequences.

#. **min_confidence**:
    * Description: specifies minimum confidence score of sequence assignment for an MA score to be calculated for the sequence.
    * Type: `float`
    * Options: any valid float in range 0-100.
    * Default: 70 (**not required**).

#. **consensus**:
    * Description: specifies whether to calculate MA from consensus spectra or individual spectra.
    * Type: `bool`
    * Options: true or false.
    * Default: true (**not required**).

#. **combine_precursors**:
    * Description: specifies whether to combine MS2 spectra from multiple unique sequence precursors before calcuating MA.
    * Type: `bool`
    * Options: true or false.
    * Default: false.

   .. note::
      * If true, this can only be done via consensus spectra and therefore **consensus** must also == true.

#. **min_peak_identity**:
    * Description: specifies minimum % of spectra a peak must be found in to be included in consensus spectra.
    * Type: `float`
    * Options: any valid float in range 0-1.
    * Default: 0.7 (**not required**).

   .. note::
      * Be careful setting this value too high if **combine_precursors** == true (peak identity is likely to be lower for spectra with different precursors).

#. **ppm_window**:
    * Description: specifies window (in ppm, parts per million) for grouping precursors and MS2 product ions in consensus spectra.
    * Type: `float`
    * Options: any valid float >= 0.
    * Default: 5 (**not required**).