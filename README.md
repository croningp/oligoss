# Oligomer Soup Sequencing (OLIGOSS)

OLIGOSS is a package for _de novo_ sequencing of linear oligomers from
tandem mass spectrometry data. 

__Note__: proper use of this package requires a good working knowledge of both
the chemistry and analytical methods used to obtain results. It is __not__ a
magic bullet to solve all of your oligomer sequencing needs! Please ensure that
all raw mass spectrometry data input to the package is of good quality. All
results should be checked and validated wherever possible, especially if running
data from a new instrument and / or oligomer class.

## Installation

OLIGOSS is available through Pip (Python Package Index):

```
pip install oligoss --user

```

## Source Code

Source code can be viewedon github: https://github.com/croningp/oligoss.git

## System Requirements

OLIGOSS was developed and tested on Ubuntu 19.10, and should therefore be compatible
with any Unix OS. As of version 0.0.1, OLIGOSS is incompabitible with Windows. Windows
compatibility will be introduced in a later version.

## Dependencies

- Python (version 3.6.0 or later)

This package was written in Python 3.7, but should be compatible with version
3.6 or later.

- mzML

Mass spectrometry data must be in .mzML file format. mzML files can be
generated from a variety of vendors using Proteowizard MS Convert, which is
freely available here: http://proteowizard.sourceforge.net/download.html

-mzmlripper

Graham Keenan's mzmlripper is required for converting mzML files into JSON
format. For full documentation: https://pypi.org/project/mzmlripper/

To install mzML ripper:

```
pip install mzmlripper --user

```
## Run

To run an OLIGOSS sequencing workflow, run the following command:

```
python -m oligoss -i input_params.json -r ripper_folder -o output_folder

```
- __input_params.json__ = input parameters file

This should contain all relevant input parameters for executing a OLIGOSS
sequencing workflow (see __Input Parameters__, below).

- __ripper_folder__ = data directory. __NOTE__: this argument can either be passed in via the command line directly (as above) or specified in the input parameters file using the _data_folder_ parameter.

This folder should contain input MS data in either mzML or ripper JSON format.

- __out_dir__ = output directory. __NOTE__: this argument can either be passed in via the command line directly (as above) or specified in the input parameters file using the _output_folder_ parameter.


All output data will be dumped to this folder.


## Sequencing Workflows

There is currently only one sequencing workflow available in OLIGOSS. The
__exhaustive screening__ workflow is to be used for sequencing oligomers with
well-characterised fragmentation pathways and known monomer libraries. For
oligomer classes with poorly characterised fragmentation pathways and / or 
data with unknown monomers, knew workflows based on previous __mass difference__
screen and Kendrick Analysis will be coming shortly.

## Configuration Files and Input Parameters

The ionization and fragmentation pathway of a set of oligomers, and the spectra
observed for MS1 and MS2 hits, varies with:

1. __oligomer class__: each oligomer class has its own set of possible
ionization and fragmentation pathways. There are usually many different ways an
oligomer can fragment.

2. __instrumentation__: the ionization source and fragmentation method used
determine which fragmentation pathways dominate for an oligomer class. The
resolution, sensitivity and other scan parameters affect how many MS1 precursors
and MS2 product ions will be detected and assigned to the correct sequence.

3. __analyte properties__: the analyte matrix (e.g. pH, salts present), as well
as the monomers used in each experiment will affect not only what oligomers are
present but how they may ionize and fragment.

Fragmentation pathways are stored in polymer-specific configuration files for
each class of oligomers (see __Polymer Configs__ section). Default operating
conditions for mass spectrometers are stored in instrument configuration files
(see __Instrument Configs__ section). All pre-configured instrument settings
are used only as defaults, and can be overwritten in the input parameters file.

All parameters which vary between individual experiments are specified in the
input parameters file (see __Input Parameters__ section, below).

<img src="img/typical_OLIGOSS_execution.png" width = 350 height = 350>

### Input Parameters

Input parameter files must be in JSON format. The following parameters are
required upon execution of every sequencing experiment:

1. __mode__:
   - Description: specifies the charge of ions.
   - Type: `str`
   - Options: either "pos" or "neg" for positive and negative ion mode,
   respectively.
2. __monomers__:
   - Description: monomers used in experiment.
   - Type: `List[str]`
   - Options: any combination of monomer one letter codes from polymer-specific configuration files.
3. __screening method__:
   - Description: specifies which sequencing workflow to use.
   - Type: `str`
   - Options: currently only option is "exhaustive".
4. __polymer_class__:
   - Description: specifies the backbone of oligomers being fragmented.
   - Type: `str`
   - Options: any valid alias associated with a polymer config file (see
   __Polymer Configs__).
5. __silico__:
   - Description: parameters that define _in silico_ properties.
   - Type: `dict`
   - Options: many sub-parameters must be defined (see __Silico Parameters__).
6. __extractors__:
   - Description: parameters that define properties essential to matching and
   filtering in MS/MS data.
   - Type: `dict`
   - Options: many sub-parameters must be defined (see __Extractor Parameters__).
7. __postprocess__:
   - Description: parameters that define properties relevant to assigning
   confidence scores to sequences based on extracted data.
   - Type: `dict`
   - Options: many sub-parameters must be defined (see __Postprocess Parameters__).

In addition to these parameters, optional parameters for mass spec model and
chromatography can be defined. __NOTE__: if instruments are not defined,
additional parameters must be defined in __silico__, __extractors__ and
__postprocess__ (see __Instrument Configs__).

9. __instrument__:
  - Description: specifies the mass spec model used to acquire MS/MS data.
  - Type: `str`
  - Options: any valid alias for a pre-configured instrument file.
10. __chromatography__:
  - Description: specifies the chromatography used to separate products
  prior to detection via MS.
  - Type: `str`
  - Options: any valid alias for a pre-configured chromatography file.

#### Silico Parameters

Silico parameters define the properties required to construst theoretical MS1 precursors
and MS2 product ions. There are general silico parameters, as well as parameters specific
to MS1 and MS2 ions, defined in silico.ms1 and silico.ms2 subparamters, respectively.

11. __min_length__:
  - Description: specifies minimum length (in monomer units) of oligomers in
  product mixtures.
  - Type: `int`
  - Options: Any integer value. __NOTE__: length of target oligomers affects size of the sequence space to screen. Limits on sequence space for screening will depend
  on oligomer class, types of monomers, length of oligomers and computing
  resources avaialable.
  - Default: __no default available__
12. __max_length__:
  - Description: specifies maximum length (in monomer units) of oligomers in
  product mixtures.
  - Type: `int`
  - Options:Any integer value. __NOTE__: length of target oligomers affects size of the sequence space to screen. Limits on sequence space for screening will depend
  on oligomer class, types of monomers, length of oligomers and computing
  resources avaialable.
  - Default: __no default available__
13. __isomeric_targets__:
  - Description: specifies isomeric sequence pool for targetting. If specified,
  only sequences isomeric to one or more isomeric targets will be included in
  the screen.
  - Type: `List[str]`
  - Options: list of any sequence strings that would be generated from
  input monomers and length distribution.
  - Default: None (__not required__).
14. __modifications__:
  - Description: specifies targets for any covalent modifications.
  - Type: `Dict[str or int, List[str]]`
  - Options: keys specify modification targets (terminal and sidechain).
    - Terminal Keys: -1 or "-1" / 0 or "0" for terminus -1 and 0, respectively.
    - Sidechain Keys: keys must be one-letter codes for monomers with compatible
    sidechains (see __Polymer Configs__).
    - Modification Values: list of modification strings. Strings must correspond
    to valid modification alias (see __Polymer Configs__).
    - Example: {"K": ["Ole", "Pal"], 0: ["Ole"]}
  - Default: None (__not required__).
15. __ms1__:
  - Description: specifies parameters for MS1 silico ions.
  - Type `dict`
  - Options: many sub-parameters (see __MS1 Silico Parameters__).
16. __ms2__:
  - Description: specifies parameters for MS2 silico ions.
  - Type `dict`
  - Options: many sub-parameters (see __MS2 Silico Parameters__).

##### MS1 Silico Parameters

MS1 silico parameters define properties required for generating theoretical MS1 precursor ions.

17. __min_z__:
  - Description: specifies minimum __absolute__ charge of MS1 precursor ions.
  - Type: `int`
  - Options: Any integer value.
  - Default: default value specified by instrument. If no instrument is specified,
  default value == __1__.
18. __max_z__:
  - Description: specifies maximum __absolute__ charge of MS1 precursor ions.
  - Type: `int`
  - Options: Any integer value >= __min_z__.
  - Default: None (__not required__). Instrument defaults can be specified in
  instrument config file.
19. __universal_sidechain_modifications__:
  - Description: specifies whether all sidechains targeted for modification
  will be modified by one or more modifying agent.
  - Type: `bool`
  - Options: true or false.
  - Default: Default can be specified in polymer config file (see __Polymer
  Configs__). If not default is specified in config, default value == __True__.
20. __universal_terminal_modifications__:
  - Description: specifies whether all termini targeted for modification
  will be modified by one or more modifying agent.
  - Type: `bool`
  - Options: true or false.
  - Default: Default can be specified in polymer config file (see __Polymer
  Configs__). If not default is specified in config, default value == __True__.
21. __max_neutral_losses__:
  - Description: specifies cap on maximum number of neutral loss fragmentation events at MS1.
  - Type: `int`
  - Options: any integer value.
  - Default: Default can be specified in instrument config file (see __Instrument
  Configs__). If default is not specified in config, default value == __None__.
22. __adducts__:
  - Description: specifies extrinsic ions present in analyte matrix that will
  affect MS1 ionization.
  - Type: `List[str]`
  - Options: list of ion strings. __NOTE__: all ions must be defined in 
  __Global Chemical Constants__.
  - Default: varies with instrument and polymer class. If no instrument- or polymer-specific defaults, this must be specified.

##### MS2 Silico Parameters

MS2 silico parameters define properties required for generating theoretical MS2 product ions.

23. __fragment series__:
  - Description: linear fragment types to be included in MS2 silico libraries.
  - Type: `List[str]`
  - Options: list of any valid MS2 fragment one letter codes in polymer config
  file (see __Polymer Configs__).
  -  Defaults: defaults vary depending on specific instrument and polymer class.
  If no default is specified in instrument config, this __MUST__ be supplied in input parameters.
24. __max_neutral_losses__:
  - Description: specifies cap on maximum number of neutral loss fragmentation
  events at MS2.
  - Type: `int`
  - Options: any integer value.
  - Default: Default can be specified in instrument config file (see __Instrument Configs__). If default is not specified in config default value == __None__.
25. __signatures__:
  - Description: specifies signature ion types to be included in MS2 silico
  libraries.
  - Type: `List[str]`
  - Options: list of any valid signature ion types in polymer config file.
  - Defaults: defaults vary depending on instrument and polymer class. If no
  default is specified in config files, default value == __None__.
26. __min_z__:
  - Description: specifies minimum __absolute__ charge of MS2 product ions.
  - Type: `int`
  - Options: Any integer value.
  - Default: default value specified by instrument. If no instrument is specified,
  no fallback default value is available so this must bepisode

Extractor parameters define properties required for screening observed MS and MS/MS data for
theoretical silico ions, and also for filtering MS data prior to screening.


28. __error__:
  - Description: specifies error threshold for matching theoretical ion _m/z_
  values to observed ions in MS/MS data.
  - Type: `float`
  - Options: Any float value between 0 and 1.
  - Default: Default value specified by instrument config. __NOTE__: no fallback
  default is available - must be specified in input parameters or instrument
  config.
29. __error_units__:
  - Description: specifies units of error tolerance value.
  - Type: `str`
  - Options: "ppm" for relative error tolerance in parts per million, or "abs"
  for absolute error tolerance in mass units (u).
  - Default: Default value specified by instrument config. __NOTE__: no fallback
  default is available - must be specified in input parameters or instrument
  config.
30. __min_rt__:
  - Description: specifies minimum retention time (in minutes) for data to
  be used for screening.
  - Type: `float`
  - Options: any valid float between 0 and length of acquisition run time.
  - Defaults: default value specified by chromatography config. If no chromatography
  is specified, default value == None (__not required__).
31. __min_rt__:
  - Description: specifies maximum retention time (in minutes) for data to
  be used for screening.
  - Type: `float`
  - Options: any valid float between 0 and length of acquisition run time.
  - Defaults: default value specified by chromatography config. If no chromatography
  is specified, default value == None (__not required__).
32. __min_ms1_total_intensity__:
  - Description: specifies minimum total intensity (taken as the sum of intensities in
  an EIC) for MS1 precursor ions to be confirmed.
  - Type: `float`
  - Options: Any float.
  - Defaults: Default value is specified by instrument config. If no instrument
  is specified, default value == None (__not required__).
33. __min_ms2_total_intensity__:
  - Description: specifies minimum total intensity (taken as the sum of intensities in
  an EIC) for MS2 product ions to be confirmed.
  - Type: `float`
  - Options: Any float >= 0.
  - Defaults: Default value is specified by instrument config. If no instrument
  is specified, default value == None (__not required__).
34. __min_ms1_max_intensity__:
  - Description: specifies minimum peak in intensity (taken as the most intensity
  signal in an EIC) for MS1 precursor ions to be confirmed.
  - Type: float
  - Options: Any float >= 0.
  - Default: Default value is specified by instrument config. If no instrument
  is specified, default value == None (__not required__).
34. __min_ms2_max_intensity__:
  - Description: specifies minimum peak in intensity (taken as the most intensity
  signal in an EIC) for MS2 product ions to be confirmed.
  - Type: float
  - Options: Any float >= 0.
  - Default: Default value is specified by instrument config. If no instrument
  is specified, default value == None (__not required__).
35. __rt_units__:
  - Description: retention time units in mzML files.
  - Type: `str`
  - Options: either "min" or "sec" for minutes and seconds, respectively.
  - Defaults: This value depends on the mass spec manufacturer, and should
  generally only be used in instrument config files. __WARNING__: please do
  not define this value in input parameters if possible. If you are unsure about the
  retention time units of your mass spec vendor, convert your mzML files to mzmlripper
  format directly and check the recorded retention times of the first and final
  spectra.
36. __pre_screen_filters__:
  - Description: specifies retention time and intensity thresholds for pre-filtering
  ripper data before screening.
  - Type: `Dict[str, float]`
  - Options: several subparameters to define (see __Pre-Screen Filters__, below).
  - Default: None (__not required__).

##### Pre-Screen Filters

Pre-screen filters are used to remove irrelevant data from raw spectra before screening.

37. __min_ms1_max_intensity__:
  - Description: specifies minimum peak in intensity in raw MS1 spectra for spectra to
  be included in later screening.
  - Type: `float`
  - Options: any valid float >= 0.
  - Defaults: None (__not required__).
38. __min_ms2_max_intensity__:
  - Description: specifies minimum peak in intensity in raw MS2 spectra for spectra to
  be included in later screening.
  - Type: `float`
  - Options: any valid float >= 0.
  - Default: None (__not required__).
39. __min_ms1_total_intensity__:
  - Description: specifies minimum total intensity in raw MS1 spectra for spectra to be
  included in later screening.
  - Type: `float`
  - Options: any valid float >= 0.
  - Default: None (__not required__).
40. __min_ms2_total_intensity__:
  - Description: specifies minimum total intensity in raw MS2 spectra for spectra to be
  included in later screening.
  - Type: `float`
  - Options: any valid float >= 0.
  - Default: None (__not required__).
41. __min_rt__:
  - Description: specifies minimum retention time for spectra to be included in later
  screening.
  - Type: `float`
  - Options: any valid float >= 0.
  - Default: None (__not required__).
42. __max_rt__:
  - Description: specifies maximum retention time for spectra to be included in later
  screening.
  - Type: `float`
  - Options: any valid float >= 0.
  - Default: None (__not required__).
43. __min_ms2_peak_abundance__:
  - Description: specifies the minimum relative abundance for most intense sequence in an MS2 spectrum. If the relative intensity of the most intense match is less than __min_ms2_peak_abundance__, any fragments found in this spectrum will be discarded.
  - Type: `float`
  - Options: any float in range 0-100.

#### Postprocess Parameters

Postprocess parameters are used to define properties relevant for assigning sequence confidence
based on observed data and how it compares to silico data, and also for any postprocessing
such as molecular assembly calculations and spectrum plots.

44. __core_linear_series__:
  - Description: core linear fragment types that are to be used to assign confidence
  for confirmed sequences from MS2 data.
  - Type: `List[str]`
  - Options: list of any valid fragment types in polymer config file. __NOTE__:
  all fragment types in __core_linear_series__ must also be present in __silico.ms2.fragment_series__.
  However, not all fragments in __silico.ms2.fragment_series__ have to be in __core_linear_series__.
  Fragment types specified _in silico_ but not in __core_linear_series__ will still be screened
  and recorded, but not used to assign confidence.
45. __exclude_fragments__:
  - Description: specifies fragments that are to be excluded from confidence calculations,
  even if they have been confirmed.
  - Type: `List[str]`
  - Options: List of any valid fragment ids in __core_linear_series__ (see below).
  - Defaults: None (__not required__).
46. __optional_core_fragments__:
  - Description: core fragments that are to be used in confidence calculations _only_
  if confirmed. If absent from extracted data, these fragments will be ignored in
  confidence assignments and final confidence score will not be lower as a result of
  their absence. However, if present, they will be used to calculate confidence.
  - Type: `List[str]`
  - Options: list of any fragment ids from core fragments.
  - Default: Default depends on instrument and polymer class. If no default is
  specified in instrument config, default value == None (__not required__).
47. __essential_fragments__:
  - Description: fragments that are essential for assigning a sequence's confidence
  score > 0.
  - Type: `List[str]`
  - Options: list of any valid fragment ids from core linear fragments.
  - Default: default depends on instrument and polymer config. If not specified in
  configs, default value == None (__not required__).
48. __dominant_signature_cap__:
  - Description: specifies maximum confidence score for sequences missing one or more
  expected "dominant" signature ion at MS2.
  - Type: `float`
  - Options: Any float in range 0-100.
  -  Default: Default value varies by instrument and polymer class. If not specified
  in configs, default value == 0 (__not required__).
49. __subsequence_weight__:
  - Description: weighting value for mean continuous fragment coverage in final confidence
  score for confirmed sequences.
  - Type: `float`
  - Options: any float in range 0-1. 0 = confidence based entirely on % confirmed fragments,
  1 = confidence based entirely on mean continuous fragment coverage.
  - Default: No default. This must be specified in input parameters file.
50. __rt_bin__ _deprecated_:
  - Description: minimum resolution (in minutes) between peaks in MS1 EICs.
  - Type: `float`
  - Options: Any valid float between 0 and length of acquisition time. __NOTE__:
  this is no longer used. It will be useful when revisiting quantification of
  individual sequences from EICs.
51. __ms2_rt_bin__ _deprecated_:
  - Description: minimum resolution (in minutes) between peaks in MS2 EICs.
  - Type: `float`
  - Options: Any valid float between 0 and length of acquisition time. __NOTE__:
  this is no longer used. It will be useful when revisiting quantification of
  individual sequences from EICs.
52. __spectral_assignment_plots__:
  - Description: specifies whether to plot annotated MS2 spectra for confirmed sequences.#
  - Type: `bool`
  - Options: true or false. If true, plots will be saved as PNG images in output directory.
  - Default: false (__not required__).
53. __min_plot_confidence__:
  - Description: specifies minimum confidence score for a sequence's spectral assignments to
  be plotted.
  - Type: `float`
  - Options: any valid float in range 0-100.
  - Default: 70.
54. __molecular_assembly__:
  - Description: specifies parameters for calculating molecular assembly values for confirmed
  sequences.
  - Type: `dict`
  - Options: several subparameters (see __Molecular Assembly Parameters__).
  - Default: several defaults for subparameters.

##### Molecular Assembly Parameters

Defines properties relevant for calculating Molecular Assembly for confirmed sequences.

55. __min_confidence__:
  - Description: specifies minimum confidence score of sequence assignment for an MA
  score to be calculated for the sequence.
  - Type: `float`
  - Options: any valid float in range 0-100.
  - Default: 70 (__not required__).
56. __consensus__:
  - Description: specifies whether to calculate MA from consensus spectra or individual spectra.
  - Type: `bool`
  - Options: true or false.
  - Default: true (__not required__).
57. __combine_precursors__:
  - Description: specifies whether to combine MS2 spectra from multiple unique sequence precursors
  before calcuating MA.
  - Type: `bool`
  - Options: true or false. __NOTE__: if true, this can only be done via consensus spectra and
  therefore __consensus__ must also == true.
  - Default: false.
58. __min_peak_identity__:
  - Description: specifies minimum % of spectra a peak must be found in to be included in consensus
  spectra.
  - Type: `float`
  - Options: any valid float in range 0-1. __NOTE__: be careful setting this value too high
  if __combine_precursors__ == true (peak identity is likely to be lower for spectra with different precursors).
   - Default: 0.7 (__not required__).
59. __ppm_window__:
  - Description: specifies window (in ppm, parts per million) for grouping precursors and
  MS2 product ions in consensus spectra.
  - Type: `float`
  - Options: any valid float >= 0.
  - Default: 5 (__not required__).

### Polymer Configs

Polymer-specific configuration files define the full scope of possible ionization and fragmentation
for an oligomer class. A proper configuration file should contain all information required to
generate full possible sequence libraries for target oligomer class, as well as MS1 precursors
and MS2 product ions. It should be divided into three sections:

1. __General Polymer Properties__:

This defines properties required for generating MS1 precursor libraries.

2. __MS2 Fragmentation Properties__:

This defines properties required for generating MS2 product ions.

3. __Modifications__:

This defines any covalent modifications and their appropriate targets.

#### General Polymer Properties

1. __MONOMERS__:
  - Description: defines monomer one-letter codes, their associated monoisotopic neutral masses
  and reactive functional groups.
  - Type: `Dict[str, list]`
  - Options: N/A
  - Example: {"A": [89.04768, [["amine", 1], ["carboxyl", 1]], "alanine"]} defines the amino acid
  monomer alanine, with a neutral monoisotopic mass of 89.04768, 1 reactive amine and 1 reactive
  carboxylic acid.
2. __MASS_DIFF__:
  - Description: defines mass lost (or gained) upon addition of a monomer to an elongating chain.
  - Type: `float` or `str`
  - Options: any valid float corresponding to neutral monoisotopic mass difference in standard
  mass units (u) _or_ a functional group string (e.g. "H2O") corresponding to the mass difference.
3. __ELONGATION__:
  - Description: defines standard number of monomer additions per elongation event (in monomer
  units).
  - Type: `int`
  - Options: any valid integer >= 1.
4. __REACTIVITY_CLASSES__:
  - Description: defines cross-reactivity of monomer functional groups.
  - Type: `Dict[str, List[list]]`
  - Options: keys must correspond to functional groups found in monomer library.
  - Example: The following example defines cross-reactivity of the "amine" functional group for
  an oligomer class, which in this case can react with either "carboxyl" or "hydroxyl" groups. The
  "amine" functional group is found in monomers "A", "B" and "C": 
              {"amine": ["carboxyl", "hydroxyl"], ["A", "B", "C"]}
5. __SYMMETRY__:
  - Description: determines whether termini are functionally equivalent (i.e. whether forward
  sequence == reverse sequence).
  - Type: `bool`
  - Options: true or false.
  - Example: __SYMMETRY__ = true for peptides as they have distinct C- and N-termini.
6. __LOSS_PRODUCTS__:
  - Description: specifies any neutral loss fragmentation events that can occur for specific monomer
  sidechains. __NOTE__: it is assumed that these fragmentation events can either be a product of
  in-source CID at MS1 or standard MS2 CID.
  - Type: `Dict[str, list]`
  - Options: keys must correspond to valid monomer one-letter codes, values must be lists of either
  monoisotopic neutral mass losses (`float`) or functional group strings corresponding to these neutral losses.
  - Example: {"N": ["NH3", "H2O"]} for monomer "N" with possible neutral loss fragmentations corresponding
  to loss of either an ammonia ("NH3") or water ("H2O") mass.
7. __IONIZABLE_SIDECHAINS__:
  - Description: specifies non-backbone ionization sites that occur at specific monomer sidechains.
  - Type: `Dict[str, dict]`
  - Options: keys must correspond to valid monomer one-letter codes. Values must be dictionaries with
  keys "pos" and "neg" defining possible ionization events for positive and negative mode, respectively.
  - Example: In the following example, monomer "K" can be ionized via proton addition in positive mode
  and the monomer "D" can be ionized via proton abstraction in negative mode: 
            {"K": {
                "pos":
                  ["H", 1, 1],
                "neg": null
              },
              "D": {
                "pos": null,
                "neg": ["-H", 1, 1]
              }}
8. __INTRINSICALLY_CHARGED_MONOMERS__:
  - Description: this defines monomers which have an intrinsic, non-exchangeable charge.
  - Type: `Dict[str, int]`
  - Options: keys must be valid monomer one-letter codes, with values equivalent to intrinsic
  charge state due to non-exchangeable ions.
  - Example: for a monomer "Z" with intrinsic charge of -2: {"Z": -2}
9. __SIDE_CHAIN_CROSSLINKS__:
  - Description: this defines any monomer-monomer crosslinks that can occur, and their effects
  on MS1 ionization and MS2 fragmentation.
  - Type: `Dict[str, dict]`
  - Options: keys must correspond to valid monomer one-letter codes. Key-Value pairs in the subdict are as follows:
    - _monomers_:
      - Description: defines other monomers that can form sidechain crosslinks with target monomer.
      - Type: `List[str]`
      - Options: list of valid monomer one-letter codes.
    - _crosslink_massdiff_:
      - Description: defines mass lost or gained upon crosslinking.
      - Type: `float` or `str`
      - Options: either a float corresponding to neutral monoisotopic mass diff or a string representing
      functional group mass diff.
    - _permissible_crosslink_charges_:
      - Description: defines permissible charge states for crosslinked moiety at the sidechain(s) of crosslinked
      monomers.
      - Type: `List[int]`
      - Options: list of any valid integer corresponding to permissible charge states.
    - _disrupt_ms2_:
      - Description: specifies whether crosslinking event disrupts standard linear fragmentation along backbone.
      - Type: `bool`
      - Options: true or false.
  - Example: The following example defines cross-linking events for the monomer "K" which, in its
  non-crosslinked state can be ionized at its sidechain (see __IONIZABLE_SIDECHAINS__). It can crosslink
  with monomers "E" and "D" via sidechain links. However, this type of crosslinking event does not disrupt
  standard linear MS2 fragmentation pathways.
    {"K": {"monomers": ["E", "D"], "crosslink_massdiff": "H2O", "permissible_crosslink_charges": [0], "disrupt_ms2": false}}

##### MS2 Fragmentation Properties

MS2 fragmentation properties are required for defining possible MS2 fragmentation pathways for an
oligomer class. This includes both linear fragment series and signature ion fragments.

10. __FRAG_SERIES__:
  - Description: dictionary that defines all linear fragmentation pathways. Linear fragmentation
  pathways are defined as any fragment series indexed stepwise along the oligomer backbone.
  - Type: `Dict[str, dict]`
  - Options: keys must correspond to fragment series one-letter codes. Properties of individual
  fragment series are defined in subparameters (see __Defining FRAG_SERIES__, below).
  - Example: see __Defining  FRAG_SERIES__ section.
11. __MS2_SIGNATURE_IONS__:
  - Description: this defines any monomer-specific signature ions that may occur. __NOTE__: the same
  fragmentation events that produce linear fragment series can also produce signature ions. However, OLIGOSS
  considers these as separate events due to the diversity of possible signature ions.
  - Type: `Dict[str, list]`
  - Options: Keys correspond to signature ion str code, values lists of monomer one-letter codes and corresponding
  free signature m/z values.
  - Example: {"Im": ["F", 120.0813], ...} defines "Im" signature fragment for monomer "F" with m/z 120.0813.
12. __MODIFICATIONS__:
  - Description: defines any covalent modifications and possible modification sites.
  - Type: `Dict[str, dict]`
  - Options: Keys must be strings corresponding to modification three-letter codes. Values are subdicts
  defining modification properties (see __Defining MODIFICATIONS__, below).

###### Defining FRAG_SERIES

The __FRAG_SERIES__ dict is used to define properties relevant to linear fragment series (i.e. fragment series that
are indexed stepwise along the oligomer backbone).

13. __default_linear__:
  - Description: specifies default linear fragmentation pathways to be included in silico libraries depending
  on the mass spec fragmentation method used to acquire data.
  - Type: `Dict[str, List[str]]`
  - Options: keys must correspond to valid fragmentation methods defined in __Instrument_Configs__. Values are list
  of linear fragment series codes for linear fragment series that are produced via the specified fragmentation method.
  - Example: In the case of "HCD" fragmentation producing "a", "b" and "y" MS2 fragments: {"HCD": ["b", "y", "a"]}.
  - __NOTE__: there is redundancy with __Instrument Configs__. Linear fragment series can also be specified for
  individual oligomer classes in instrument config files. These can also be overwritten directly in input parameters
  file.
14. __default_core__:
  - Description: specifies default core linear fragment series (i.e. linear series used in confidence assignments)
  depending on the mass spec fragmentation method used to acquire data.
  - Type: `Dict[str, List[str]]`
  - Options: keys must correspond to valid fragmentation methods defined in __Instrument Configs__. Values are lists
  of linear fragment series codes for linear fragment series that are produced via the specified fragmentation method
  _and_ are required for assigning confidence scores.
  - Example: In the case of "HCD" fragmentation producing core fragment series "b" and "y": {"HCD": ["b", "y"]}.
  - __NOTE__: there is redundancy with __Instrument Configs__. Linear fragment series can also be specified for
  individual oligomer classes in instrument config files. These can also be overwritten directly in input parameters
  file.
15. __terminus__:
  - Description: specifies "home" terminus from which linear fragment series is indexed.
  - Type: `int`
  - Options: either 0 or -1 for fragment series indexed from terminus 0 and -1, respectively.
  - __NOTE__: for oligomer classes with __symmetry__ == False, terminus is irrelevant.
16. __mass_diff__:
  - Description: specifies _neutral_ mass difference between a fragment and its corresponding intact neutral
  sequence slice.
  - Type: `str` or `float`
  - Options: valid float corresponding to mass difference in mass units (u) or functional group string representing
  a neutral monoisotopic mass corresponding to mass (e.g. "OH", "H2O").
17. __fragmentation_unit__:
  - Description: specifies increment of fragment indices when producing linear fragment series.
  - Type: `Dict[str, [int or str]`
  - Options: must be a key for "pos", "neg" defining fragmentation unit in positive and negative mode, respectively.
  Values must either be ints >= 0 or strings representing int value (most commonly "ELONGATION_UNIT" if
  __fragmentation_unit__ == __ELONGATION_UNIT__).
  - Example: __fragmentation_unit__ will equal 1 or __ELONGATION_UNIT__ for the majority of oligomer classes. A possible exceptions to this would be for alternating copolymers with alternating backbone links:
18. __start__:
  - Description: start position of fragment series relative to __terminus__.
  - Type: `int`
  - Options any valid integer >= 0.
  - Example: 0 for a fragment series that begins immediately at __terminus__, 1, 2, 3 for fragment series that
  begins 1, 2 or 3 indices away from __terminus__.
19. __end__:
  - Description: end position of fragment series relative to other terminus (i.e. terminus 0 and -1 for __terminus__
  == 1 and __terminus__ == 0, respectively).
  - Type: `int`
  - Options: any valid integer >= 0. 
  - Example: 0 for a fragment series that terminates at final index on backbone, 1, 2, 3 for fragment series that
  terminates 1, 2 or 3 indexes away from final index on backbone.
20. __intrinsic_charge__:
  - Description: defines any _non-exchangeable_ ions associated with fragments of a particular series.
  - Type: `Dict[str, int]`
  - Options: keys must be "pos" and "neg" for positive and negative mode, respectively. Values must be
  integers representing intrinsic charge value.
  - Example: {"pos": 1, "neg": null} for a fragment series with intrinsic charge of 1 in positive mode but no
  intrinsic charge in negative mode.
  - __NOTE__: do not confuse this with __intrinsic_adduct__. By definition MS2 fragment series are charged by default. However, this can be a result of either _non-exchangeable_ or _exchangeable_ ions. __intrinsic_charge__
  defines charge state due to _non-exchangeable_ ions.
  - __NOTE__: __NOT_REQUIRED__. This property does not need to be defined if __intrinsic_adduct__ is defined. However,
  at least one of these properties must be defined to account for fragment charge.
21. __intrinsic_adduct__:
  - Description: defines any _exchangeable_ ions associated with fragments of a particular series.
  - Type: `Dict[str, str]`
  - Options: keys must be "pos" and "neg" for positive and negative mode, respectively. Values must correspond
  to strings representing adducts stored in __Global Chemical Constants__.
  - Example: {"pos": "H", "neg": "-H"} for a fragment series that is intrinsically protonated in positve mode but
  deprotonated in negative mode.
  -__NOTE__: this property should only be used to define _exchangeable ions_ (i.e. ions that can be swapped for
  extrinsic ions in sample matrix). Do not confuse with _non-exchangeable ions_, which are defined in  __intrinsic_charge__.
22. __exceptions__:
  - Descriptions: for oligomer classes with mixed backbones (i.e. more than one backbone bond type that can be
  fragmented at MS2), fragmentation properties may differ depending on what type of bond is being fragmented at
  a particular index.
  - Type: `Dict[str, dict]`
  - Options: keys must be "pos" and "neg" to define exceptions to standard fragmentation rules in positive and negative mode, respectively. Values define exceptions to any combination of previously described MS2 fragmentation properties for linear fragment series.
  - Format: {mode (str): {func_group: {prop: {"positions": List[int], "start": int, "end": int, "exception_value: Value}}}}
    - _mode_ == either "pos" or "neg" for positive or negative mode.
    - _func_group_ == functional group that causes exception to standard fragmentation pathway.
    - _prop_ == the property for which the exception may apply.
    - _positions_: defines list of indexes in subsequence at which exception applies. Some fragmentation exceptions
    only apply when the non-standard backbone link is in a particular position in the fragment subsequence.
    - _start_: defines start position at which exception applies, relative to home terminus.
    - _end_: specifies number of indices away from end terminus at which exception no longer applies
    - _exception_value_: the substituted value to use for the property if exception applies.
  - Example: The following example is for a fragment series with exception to __mass_diff__ in cases where a bond between a "hydroxyA"-containing monomer is being fragmented. The exception applies when the "hydroxyA"-containing monomer occurs at the final index of the subsequence. The exception applies from the very first index of the fragment series but ends one index away from the end terminus:
  {"pos": {"hydroxyA": {"mass_diff": {"positions": [-1], "start": 0, "end": 1, "exception_value": 26.98709}}}}

### Instrument Configs

Instrument configuration files are used to store information on mass spectrometers used routinely for experiments. These can define resolution (in terms of error tolerance for matching peaks), sensitivity (in terms of minimum intensity thresholds for detection), and fragmentation methods.

1. __error__:
  - Description: this defines the default error threshold for an instrument when matching peaks.
  - Type: `float`
  - Options: any valid float >= 0. This can correspond to relative error threshold (parts per million, ppm) or absolute error threshold (mass units, u).
2. __error_units__:
  - Descriptions: specifies units of default __error__.
  - Type: `str`
  - Options: either "ppm" or "abs" for relative and absolute error thresholding, respectively.
3. __rt_units__:
  - Description: specifies the default retention time units in mzML files. This is a vendor-specific property outside the control of OLIGOSS (e.g. mzML files generated from Bruker mass specs have retention time units of seconds, while ThermoScientific mass spec units are in minutes).
  - Type: `str`
  - Options: either "min" or "sec" for seconds or minutes, respectively.
  - __NOTE__: if you are unsure about the retention time units in your mzML files, this is usually not specified in the raw mzML itself. Spectra in output rippers are sorted by retention time, so it should be straightforward to work out __rt_units__ for your mass spec from the recorded retention times of the first and last spectra (assuming you know total acquisition time).
4. __min_ms1_max_intensity__:
  - Description: specifies default minimum peak in intensity for accepting an MS1 EIC as valid.
  - Type: `float`
  - Options: any valid float >= 0.
5. __min_ms2_max_intensity__:
  - Description: specifies default minimum peak in intensity for accepting an MS2 EIC as valid.
  - Type: `float`
  - Options: any valid float >= 0.
6. __fragmentation__:
  - Description: specifies fragmentation methods available at every stage of tandem mass spectrometry for an instrument.
  - Type: `Dict[str, List[str] or str]`
  - Options: keys must include "ms1", "ms2" and (optionally) "msn" for defining fragmentation methods at MS1, MS2 and MS3+ levels respectively.
  - Example: for a mass spec with "neutral" fragmentation (i.e. is-CID) at MS1, and "HCD" and "CID" at MS2-n: 
  {"ms1": "neutral", "ms2": ["HCD", "CID", "neutral"], "msn": ["HCD", "CID", "neutral"]}.
7. __pre_screen_filters__:
  - Description: specifies default intensity thresholds for pre-filtering spectra before screening.
  - Type: `Dict[str, float]`
  - Options: keys:
    - _min_ms1_max_intensity_: specifies minimum intensity of base peak for MS1 spectra to be included in screening.
    - _min_ms2_max_intensity_: specifies minimum intensity of dominant ion for MS2 spectra to be included in screening.
8. __polymer_classes__:
  - Description: defines default fragmentation ionization, fragmentation and some postprocessing parameters for individual oligomer classes when using the instrument.
  - Types: `Dict[str, dict]`
  - Options: keys must be valid polymer config aliases. Values are subdicts defining default parameters for __silico_ms1__, __silico_ms2__, __extractors__ and __postprocessing__.