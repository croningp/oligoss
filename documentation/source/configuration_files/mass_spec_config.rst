.. _Instrument-Config:

##############################
Instrument Configuration Files 
##############################

All parameters relevant to extraction and processing of data can be stated explicitly in 
input parameter files (see `input parameters<Input-Parameters>`). However, sensible default
values for many of these parameters can be determined based on the particular mass spec model
used. Therefore, OLIGOSS allows users to construct instrument-specific configuration files.

.. _Mass-Spec-Config:

Mass Spec Configuration
=======================

To avoid repeatedly defining the same parameters for each new sequencing experiment, 
default values for the following parameters can be saved in mass spec configuration files
for each instrument:

#. **error**:
   * Description: this defines the default error threshold for an instrument when matching peaks.
   * Type: `float`
   * Options: any valid float >= 0. This can correspond to relative error threshold (parts per million, ppm) or absolute error threshold (mass units, u).


#. **error_units**:
   * Description: specifies units of default **error**.
   * Type: `str`
   * Options: either "ppm" or "abs" for relative and absolute error thresholding, respectively.


#. **rt_units**:
   * Description: specifies the default retention time units in mzML files. This is a vendor-specific property outside the control of OLIGOSS (e.g. mzML files generated from Bruker mass specs have retention time units of seconds, while ThermoScientific mass spec units are in minutes).
   * Type: `str`
   * Options: either "min" or "sec" for minutes and seconds, respectively.
   * **NOTE**:  if you are unsure about the retention time units in your mzML files, this is usually not specified in the raw mzML itself. Spectra in output rippers are sorted by retention time, so it should be straightforward to work out rt_units for your mass spec from the recorded retention times of the first and last spectra (assuming you know total acquisition time).


#. **min_ms1_max_intensity**:
   * Description: specifies default minimum peak in intensity for accepting an MS1 EIC as valid.
   * Type: `float`
   * Options: any valid float >= 0.


#. **min_ms2_max_intensity**:
   * Description: specifies default minimum peak in intensity for accepting an MS2 EIC as valid.
   * Type: `float`
   * Options: any valid float >= 0.


#. **fragmentation**:
   * Descroption: specifies fragmentation methods available at every stage of tandem mass spectrometry for an instrument.
   * Type: `Dict[str, List[str]`]
   * Options: keys must include "ms1", "ms2" and (optionally) "msn" for defining fragmentation methods at MS1, MS2 and MS3+ levels respectively.
   * Example for an ESI-mass spec with neutral loss fragmentation (i.e. is-CID) at MS1, and "HCD" and "CID" at MS2-n::
      {"ms1": "neutral", "ms2": ["HCD", "CID", "neutral"], "msn": ["HCD", "CID", "neutral"]}


#. **pre_screen_filters**:
   * Description: any valid parameter from `pre-screen filters<Pre-Screen-Filters>`


#. **polymer_classes**:
   * Description: any silico or postprocess parameter for specific oligomer class.
