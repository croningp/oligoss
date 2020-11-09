.. _Postprocess-Parameters:

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
