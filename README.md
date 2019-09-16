# PolymerMassSpec

PolymerMassSpec is a python project designed to enable _de novo_ sequencing and
analysis of highly heterogeneous mixtures of small to medium-sized linear
polymers from mass spectrometry data.

__Note__: proper use of this package requires a decent working knowledge of
both the chemistry and analytical methods used to obtain results. It is
__not__ a magic bullet that can solve all of your complex polymer soup problems!
Please ensure that all results are validated wherever possible, and that
the raw mass spectrometry data input to the package is of good quality. See
our __"Validation Tips and Tricks"__ section below for advice on how to do
this.

## Dependencies

- Python (version 3.6.0 or later)

This package was written in Python 3.7, but should be compatible with version
3.6 or later.

- mzML

Mass spectrometry data must be in .mzML file format. mzML files can be
generated from a variety of vendors using Proteowizard MS Convert, which is
freely available here: http://proteowizard.sourceforge.net/download.html

- mzmlripper

mzML files must be converted to JSON using Graham Keenan's mzmlripper,
instructions and documentation for which can be found here:

http://datalore.chem.gla.ac.uk/Origins/mzmlripper.git

## Installation

Repo can be cloned from:

```
http://datalore.chem.gla.ac.uk/DBD/polymermassspec.git

```

## Run

__1__ Convert raw data files to mzmlripper JSON files

__2__ Ensure _Input Parameters_ JSON file is filled in with correct run
parameters (see __"Run Parameters"__ section)

__3__ Run executable script:

Experiments are run from the executable, which can be found here:

_polymersoup/executable.py_

To run the script:

```
python -m polymersoup.executable input_parameters_file.json

```
## Run Parameters

All run parameters are passed in through a single JSON file.

Input parameters are divided into three sections:

### 1. silico parameters:

__General silico parameters__:

```
__"mode"__:

- Description: this parameter specifies whether mass spec acquisition is in positive or negative mode

- Options: either _"pos"_ or _"neg"_ for negative and positive mode, respectively

```
__MS1 silico parameters__:

```
__"monomers"__:

- Description: list of monomers used in your reaction(s). These __MUST__ match monomer one letter codes supplied in polymer-specific config file (see __"Polymer Config File"__ section, below)

- Options: any combination of monomer one letter codes supplied in
polymer-specific config file

__"max_length"__:

- Description: maximum length of target sequence(s)

- Options: any integer

- Recommendations: we do not recommend attempting to sequence product pools with more than 2 million potential sequences at this stage. Anything up to and including this should be doable (if not, please let us know!!!)

__"min_length"__:

- Description: specifies minimum length of target sequence(s)

- Options: integer

- Recommended: 1 for most experiments

__"ms1_adducts"__:

- Description: list of adducts to screen for in MS1 data. Adducts are supplied as lists of strings, each adduct string __MUST__ be present in GlobalChemicalConstants.py. If you require adducts that are not in this file, please let us know.

- Options: any combination of adducts in the GlobalChemicalConstants file that match overall charge state of mode

- Recommended: ["H"] for most samples in positive mode

__"min_z"__:

- Description: specifies minimum __absolute__ charge state for MS1 precursor ions.

- Options: integer

- Recommended: 1 for standard samples. __NOTE__: if you have any monomers that are _intrinsically_ charged, do __NOT__ account for this by increasing "min_z". Rather, update the polymer-specific config file
and this will be factored in to any runs with charged monomers (see __"Polymer Config File"__ section, below)

__"max_z"__:

- Description: specifies maximum __absolute__ charge state for MS1 precursor ions.

- Options: integer

__"losses"__:

- Description: specifies whether to include side chain-specific neutral loss products for MS1 precursor ions.

- Options: true or false

- Recommended: true, unless you are certain from previous observations that these loss products will not be observed. __NOTE__: this can drastically affect the results, __be careful__

__"max_neutral_losses"__:

- Description: specifies whether to place an upper limit on the number of side chain-specific neutral loss products to include for any given MS1 precursor

- Options: integer or null

- Recommended: null, unless you have a problem with false positives due to screening for too many precursor ions (this may be an issue for sequence pools with a large excess of loss product-prone side chains)

__"loss_products_adducts"__:

- Description: specifies whether to add extra adducts to loss products

- Options: any list of adducts found in GlobalChemicalConstants.py

- Recommended: null

__"chain_terminators"__:

- Description: list of monomers in your reaction that terminate chain elongation for reasons __other than__ functional group cross-reactivities that have already been specified in polymer-specific config file

- Options: list of monomer one letter codes. All one letter codes must be found in polymer-specific config file

- Recommended: null, unless you know for certain that certain monomers will terminate your sequences for reasons other than their inherent reactivity. This will be used in cases where sequences have been treated with side chain-specific cleaving agents (e.g. chymotrypsin for peptides).

__"universal_rxn"__:

- Description: this specifies whether all input monomers are universally cross-reactive. If this is the case, the code only needs to perform very simple operations to generate _in silico_ data. Otherwise, more complex (and therefore more _time consuming_) operations are required to generate _in silico_ data

- Options: true or false

- Recommended: true if this true, otherwise you __MUST__ make sure that monomer reactivity classes are defined correctly in the polymer-specific config file

__"terminal_tags"__:

- Description: defines list of tags that are present at peptide termini. tags at "0" are for start terminus, "-1" are for end terminus

- Options: null or list of tag codes (__currently set up for monomer tags__, but other tags are coming next week)

__"side_chain_tags"__:

- Description: dictionary where keys are amino acids and their values consist of a list of possible side
chain modifications.

- Options: null or a dictionary of the format: key = amino acid, value = list of side chain modifications for that amino acid.

__"cyclic_sequences"__:

- Description: states whether or not cyclic sequences may be formed during the experiment.

- Options: True or False

- Recommended: False unless you are certain that cyclic sequences can be formed during your experiments.

__"isobaric_targets"__:

- Description: to be completed

```
__MS2 silico parameters__:

```
__"fragment series"__:

- Description: list of fragment series one letter codes

- Options: any list of fragment one letter codes found in polymer-specific config file

- Recommended: ["b", "y"] for peptides fragmented via CID

__"ms2_adducts"__:

- Description: list of adducts to add to MS2 fragments other than intrinsic adducts specified in polymer-specific config file

- Options: any list of adducts found in GlobalChemicalConstants

- Recommended: ["Na"] for standard experiments in positive mode

__"ms2_losses"__:

- Description: specifies whether to include side chain-specific neutral loss products for MS2 fragment ions

- Options: true or false

- Recommended: true

__"ms2_max_neutral_losses"__:

- Description: specify upper cap on number of side chain-specific neutral loss products to include for MS2 fragments

- Options: null or integer

- Recommended: null

__"ms2_loss_products_adducts"__:

- Description: specify whether to add any adducts on to MS2 neutral loss products

- Options: null or list of any adducts found in GlobalChemicalConstants.py

- Recommended: null

__"add_signatures"__:

- Description: specify whether to add monomer-specific signature ions for MS2 fragments

- Options: true or false

- Recommended: true

__"signatures"__:

- Description: if monomer MS2 signatures are to be added, this specifies any particular subsets of signatures that are to be added. If "add_signatures"==true and null is supplied, ALL signature types will be added

- Options: null or list of signature types specified in polymer-specific config file

- Recommended: null

__"min_z"__:

- Description: specifies minimum charge of MS2 ions

- Options: integer

- Recommended: 1

__"max_z"__:

- Description: specifies maximum charge of MS2 ions

- Options: integer

- Recommended: 1 for standard CID of small to medium-sized polymers, with the exception of polymers with excess of intrinsically charged monomers and / or large multi- metal centre transition metal complexes

```

## 2. extractor_parameters:

__general extractor parameters__:

```
__"error"__:

- Description: specifies error tolerance threshold for matching target ions to ions observed in mass spectra. This can be supplied as an absolute value (amu) or relative error (i.e. ppm)

- Options: float

- Recommended: 0.01 (absolute)

__"err_abs"__:

- Description: specifies whether error units are in absolute mass units or ppm

- Options: true or false

- Recommended: true

__"min_ms2_peak_abundance"__:

- Description: specifies minimum relative abundance of most intense matching MS2 peak for a target sequence in observed MS2 spectra, expressed as a % of the most intense peak observed in the spectrum

- Options: float (range = 0 to 100)

- Recommended: 90-100

__"pre_run_filter"__:

- Description: specifies whether to filter spectra before beginning sequence screening. If true, spectra will be filtered through pre_screen_filters and sequences will be screened against truncated data set of spectra that pass these filters.

- Options: true or false

- Recommended: true

```

__pre_screen_filters__:

```

__"min_rt"__:

- Description: specifies minumum retention time of raw spectra. Any spectra with retention time lower than this will be discarded from consideration

- Options: float

- Recommneded: 0 unless you know your chromatography well enough

__"max_rt"__:

- Description: specifies maximum retention time of raw spectra. Any spectra with retention time higher than this will be discarded from consideration

- Options: float

- Recommended: run time (in minutes) of your method unless you know exactly what you are doing with this

__"essential_signatures"__:

- Description: list of monomers that MUST have signatures present in MS2 spectra. Any MS2 spectra without these signatures will be discarded.

- Options: null or list of monomer one-letter codes. These monomers MUST have signature ions in polymer-specific config file.

- Recommended: null, unless you know for certain every sequence has these monomers AND  that these monomer signatures will always show up in MS2 spectra.

__"signature_types"__:

- Description: specifies list of signature codes to include in essential signatures. These MUST correspond to signature types in polymer-specific config file

- Options: any list of signature type found in polymer-specific config file

- Recommended: ["Im"] for peptides using depsipeptide config file, as this signature class represents immonium ions

__"signature_ms_level"__:

- Description: specifies which MS level signatures will be present

- Options integer

- Recommended: 2, as presently we only have data for MS2 signatures

__"massdiff_bins"__:

- Description: *to be completed*

- Options: true or false

- Recommended: *to be completed*

__"ms2_precursors"__:

- Description: list of precursor ions for matching to MS2 parents in observed MS2 spectra. Any MS2 spectra that do not have a parent matching one or more of these precursors will be discarded from consideration.

- Options: null or list of floats

- Recommended: null unless doing very targeted screening for specific small subset of products (which general screen should pick up anyway...)

__"min_MS1_total_intensity"__:

- Description: specifies minimum total intensity of MS1 spectra. Any MS1 spectra with total intensity lower than this will be discarded from consideration

- Options: null or float

- Recommended: null

__"min_MS2_total_intensity"__:

- Description: specifies minimum total intensity of MS2 spectra. Any MS2 spectra with total intensity lower than this will be discarded from consideration

- Options: null or float

- Recommended: null

__"min_MS1_max_intensity"__:

- Description: specifies minimum maximum intensity of MS1 spectra. Any MS1 spectra with maximum intensity lower than this will be discarded from consideration

- Options: null or float

- Recommended: null

__"min_MS2_max_intensity"__:

- Description: specifies minimum maximum intensity of MS2 spectra. Any MS2 spectra with maximum intensity lower than this will be discarded from consideration

- Options: null or float

- Recommended: null

```
## 3. postprocess parameters:

```

__"exclude_frags"__:

- Description: list of specific fragment ids to exclude from consideration in confidence assignments ONLY (will still be included in retention time assignment)

- Options: null or any list of valid fragment ids

- Recommended: null

__"optional_core_frags"__:

- Description: list of specific fragment ids to exclude from consideration in confidence assignments IF they are absent

- Options null or any list of valid fragment ids

- Recommended: ["b1"] for peptides in positive mode, as this particular MS2 fragment is notorious for not showing up in MS2 spectra of peptides

__"core_linear_series"__:

- Description: list of fragment series one letter codes for fragment series that should be used in confidence calculations. NOTE: it is assumed that each series in this list will receive an equal weighting in confidence calculation; if this is not suitable for your chemistry, please let us know!!!

- Options: list of valid fragment one letter codes. NOTE: this MUST not include signatures

- Recommended: ["b", "y"] for peptides in positive mode (["y"] for peptides in negative mode)

__"excluded_fragments"__:

- Description: list of specific fragment ids to exclude from confidence assignments AND retention time assignments

- Options: null or any list of specific fragment ids

- Recommended: null

__"dominant_signature_cap"__:

- Description: specifies an upper cap on confidence assignments for sequences that are expected to have abundant monomer-specific MS2 signature(s) but have one or more of these missing

- Options: float (range = 0 to 100)

- Recommended: 70

__"essential_fragments"__:

- Description: specifies specific fragment ids that MUST be present for a sequence to be assigned with confidence

- Options: null or list of specific valid fragment ids

- Recommended: null

__"subsequence_weight"__:

- Description: specifies weighting to place on continuous fragment coverage in calculating overall confidence of a sequence assignment. This is supplied as a decimal fraction.

- Options: float (range = 0 to 100)

- Recommended: we are still trying to work this out. Probably between 0.25 and 0.5

__"min_rt"__:

- Description: mininimum retention time to assign to confirmed sequences

- Options: null or float

- Recommended: null

__"max_rt"__:

- Description: maximum retention time to assign confirmed sequences

- Options: null or float

- Recommended: null

__"Rt_bin"__:

- Description: specifies minimum separation of MS1 peaks in extracted ion chromatograms (units = minutes)

- Options: float

- Recommended: 0.25 for LC-MS, 0 without chromatography

__"backup_Rt_bin"__:

- Description: specifies a back-up retention time bin if peaks cannot be assigned with Rt_bin. This may be useful in cases where particular sequence subsets do not separate well using chosen chromatography method (which happens sometimes for complex mixtures)

- Options: float (should be lower than Rt_bin)

- Recommended: 0.1

__"ms2_Rt_bin"__:

- Description: specifies maximum discrepancy between retention time of unique MS2 fragment and MS1 precursor

- Options: float

- Recommended: 0.5 (this may be too generous)

__"ms2_Rt_flexible"__:

- Description: specifies whether to allow for the possibility of a mismatch between observed MS2 retention time and MS1 precursor within specified ms2_Rt_bin. If true, ms2_Rt_bin tolerance will be widened until a retention time can be assigned.

- Options: true or false

- Recommended: we need to test this urgently

__"min_viable_confidence"__:

- Description: specifies minimum confidence score for sequences to be assigned final retention time and intensity values. Sequences with scores lower than this will not be processed for retention time and / or intensities

- Options: float

- Recommended: 55-60

__"min_relative_intensity"__:

- Description: specifies minimum relative intensity (as % intensity of most intense peak) in EICs. Any peaks with relative intensities lower than this will be discarded, which may lead to less abundant sequences not being assigned.

- Options: null or float (range = 0 to 100)

- Recommended: null, unless you are deliberately attempting to assign retention time of very low confidence sequences (which probably won't work anyway..)

__"plot EICs"__:

- Description: option to plot all EICs and save the plots as png files into your output folder.
    NOTE: not finished yet.

- Options: true or false

```
## directories:

```
__"ripper_folder"__:

- Description: specifies location of folder containing mzml ripper JSON files

- Options: string file path to folder

- Recommended:

__"output_folder"__:

- Description: specifies location of output folder, where output data will be saved

- Options: full string file path to output folder

- Recommended: scapa data folder (obviously)

```

__"Polymer Config File"__:

# SECTION 1: GENERAL MS1 AND CHEMISTRY RULES:

```
__"MONOMERS"__:

- Description: dictionary of monomer one letter codes, along with associated neutral monoisotopic masses and reactivity classes (functional groups). This dictionary should specify the neutral monoisotopic mass for each monomer, its reactive functional groups, and the number of each reactive functional group. This information is essential to ensure that all theoretical sequences generated are chemically feasible.

- Format:  {'X': [mass, [[rxn_class1, n], [rxn_class2, y]]}
where X = monomer one letter code; mass = neutral monoisotopic mass; rxn_class1 and rxn_class2 = reactivity classes (functional groups - e.g. amine, aldehyde); n = number of rxn_class1 functional groups; y = number of rxn_class2 functional groups.

__"MASS_DIFF"__:

- Description: the mass difference when adding an additional monomer on to a polymer chain.

- Format: float

- Options: For condensation polymer = H2O

__"ELONGATION_UNIT"__:

- Description: the number of additional monomer units typically added when elongating a polymer.

- Format: integer

- Options: This will typically be 1 for most polymers.

__"REACTIVITY_CLASSES"__:

- Description: dictionary of reactivity classes with associated compatible classes and monomers.

- Format: {'classA' : [['classX', 'classY'], ['A', 'B']}
where classA = reaction class, classX and classY = classes that are cross-reactive with classA, and A and B are monomers within reaction class classA.

__"SYMMETRY"__:

- Description: bool to define whether polymer is identical at both ends.

- Options: true or false

- Recommended: Set to false for polymers with different termini (e.g. N- and C- termini for peptides), true for polymers with identical functional groups at both termini in linear chains.

__"CHAIN_TERMINATORS"__:

- Description: list of monomers that terminate chain elongation.

- Options: null unless you are certain a certain monomer will stop elongation of the chain in your reactions.

__"LOSS_PRODUCTS"__:

- Description: dictionary of monomer one letter codes and associated side chain neutral loss products, i.e. masses that can be lost from the monomer side chain.

__"IONIZABLE_SIDECHAINS"__:

- Description: dictionary of monomers that can be ionised with extra adducts at the SIDE CHAIN, with associated adducts, minimum and maximum absolute charge states in both positive and negative mode.

- Format: {"X": {"pos": (adduct, a, b), "neg": (adduct, a, b)}}
where "X" = monomer one letter code, adduct = adduct string (must be found in either CATIONS OR ANIONS in GlobalChemicalConstantss), a = min side chain charge, b = max side chain charge for ionized form.

__"INTRINSICALLY_CHARGED_MONOMERS"__:

- Description: dictionary of monomers that have an intrinsic charge (i.e. charged without addition of adducts), with associated lists of permissible adducts

```

# SECTION 2: MS2 FRAGMENTATION

- Description: This section should contain all the information required to construct a basic MS2 fragment series for linear polymers. Each fragment type is defined as a key in the FRAG_SERIES dict. Typically, fragment keys are one letter codes used to denote fragments - e.g. the standard 'b' and 'y' fragment series for peptides.
- When building fragment series, the fragment generator will use the following convention: f'{frag}{n}'
frag = fragment one letter code, n = number of monomers in fragment (e.g. b1, b2, b3, y1, y2, y3 etc...)


## __"FRAG_SERIES"__:

- Description: dictionary of fragment series one letter codes and associated properties.

- Format: Each fragment series key = fragment one letter code, value = subdictionary containing all relevant fragment properties.

```
__"terminus"__:

- Description: end in which fragmentation starts.

- Options: -1 (for end of the sequence) or 0 (for start of the sequence)

__"mass_diff"__:

- Description: the difference in mass between the fragment and corresponding subsequence. This may vary depending on whether cations or anions are being fragmented, therefore separate mass_diffs are specified for positive mode ('pos') and negative mode ('neg') mass spec. Example: the 'y3' fragment of peptide sequence 'AGVS' = the mass of (GVS+H) in positive mode, and (GVS-H) in negative mode.

__"fragmentation_unit"__: 

- Description: the minimum number of monomer units typically added and / or removed at a time when building a fragment series.

- Recommended: Default is ELONGATION_UNIT.

```

## __"MS2_SIGNATURE_IONS"__:

```
- Description: MS2 fragments which can be used as markers for monomers and / or small subsequences.

- Format: nested dictionary of format {"type": {{"X"}: [signature_mass]}}
where "type" is the signature ion type, "X" is the monomer and signature_mass is a float which corresponds to the mass of the signature ion for the monomer ("X").




