.. _Extractor-Parameters:

Extractor Parameters
====================

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

.. _Pre-Screen-Filters:

Pre-Screen Filters
------------------

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
