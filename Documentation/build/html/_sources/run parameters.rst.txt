Run Parameters
##############

All run parameters are passed in through a single JSON file.

Input parameters are divided into three sections:

1. silico parameters
2. extractor parameters
3. post processing parameters

Silico Parameters
=================

general silico parameters
-------------------------
- "mode":

    - **Description:** this parameter specifies whether mass spec acquisition is in positive or negative mode.

    - **Options:** either "pos" or "neg" for negative and positive mode, respectively.


MS1 silico parameters
---------------------
- "monomers":

    - **Description:** list of monomers used in your reaction(s). These **MUST** match monomer one letter codes supplied in polymer-specific config file (see - **"Polymer Config File"** section).

    - **Options:** any combination of monomer one letter codes supplied in polymer-specific config file.

- "max_length":

    - **Description:** maximum length of target sequence(s).

    - **Options:** any integer.

    - **Recommended:** we do not recommend attempting to sequence product pools with more than 2 million potential sequences at this stage. Anything up to and including this should be doable (if not, please let us know!!!).

- "min_length" :

    - **Description:** specifies minimum length of target sequence(s)

    - **Options:** integer

    - **Recommended:** 1 for most experiments

- "ms1_adducts" :

    - **Description:** list of adducts to screen for in MS1 data. Adducts are supplied as lists of strings, each adduct string **MUST** be present in GlobalChemicalConstants.py. If you require adducts that are not in this file, please let us know.

    - **Options:** any combination of adducts in the GlobalChemicalConstants file that match overall charge state of mode.

    - **Recommended:** ["H"] for most samples in positive mode.

- "min_z" :

    - **Description:** specifies minimum **absolute** charge state for MS1 precursor ions.

    - **Options:** integer.

    - **Recommended:** 1 for standard samples. **NOTE**: if you have any monomers that are _intrinsically_ charged, do **NOT** account for this by increasing "min_z". Rather, update the polymer-specific config file and this will be factored in to any runs with charged monomers (see - **"Polymer Config File"** section, below).

- "max_z" :

    - **Description:** specifies maximum **absolute** charge state for MS1 precursor ions.

    - **Options:** integer.

- "losses" :

    - **Description:** specifies whether to include side chain-specific neutral loss products for MS1 precursor ions.

    - **Options:** true or false.

    - **Recommended:** true, unless you are certain from previous observations that these loss products will not be observed. **NOTE**: this can drastically affect the results, **be careful**.

- "max_neutral_losses" :

    - **Description:** specifies whether to place an upper limit on the number of side chain-specific neutral loss products to include for any given MS1 precursor.

    - **Options:** integer or null.

    - **Recommended:** null, unless you have a problem with false positives due to screening for too many precursor ions (this may be an issue for sequence pools with a large excess of loss product-prone side chains).

- "loss_products_adducts" :

    - **Description:** specifies whether to add extra adducts to loss products.

    - **Options:** any list of adducts found in GlobalChemicalConstants.py.

    - **Recommended:** null.

- "chain_terminators" :

    - **Description:** list of monomers in your reaction that terminate chain elongation for reasons **other than** functional group cross-reactivities that have already been specified in polymer-specific config file.

    - **Options:** list of monomer one letter codes. All one letter codes must be found in polymer-specific config file.

    - **Recommended:** null, unless you know for certain that certain monomers will terminate your sequences for reasons other than their inherent reactivity. This will be used in cases where sequences have been treated with side chain-specific cleaving agents (e.g. chymotrypsin for peptides).

- "universal_rxn" :

    - **Description:** this specifies whether all input monomers are universally cross-reactive. If this is the case, the code only needs to perform very simple operations to generate *in silico* data. Otherwise, more complex (and therefore more **time consuming**) operations are required to generate *in silico* data.

    - **Options:** true or false.

    - **Recommended:** true if this true, otherwise you **MUST** make sure that monomer reactivity classes are defined correctly in the polymer-specific config file.

- "terminal_tags" :

    - **Description:** defines list of tags that are present at peptide termini. Tags at "0" are for start terminus, "-1" are for end terminus.

    - **Options:** null or list of tag codes.

- "side_chain_tags" :

    - **Description:** dictionary where keys are amino acids and their values consist of a list of possible side chain modifications.

    - **Options:** null or a dictionary of the **Format:** key = amino acid, value = list of side chain modifications for that amino acid.

- "cyclic_sequences" :

    - **Description:** states whether or not cyclic sequences may be formed during the experiment.

    - **Options:** true or false.

    - **Recommended:** False unless you are certain that cyclic sequences can be formed during your experiments.

- "isobaric_targets" :

    -  **Description:** to be completed.


MS2 silico parameters
---------------------


- "fragment series" :

    - **Description:** list of fragment series one letter codes.

    - **Options:** any list of fragment one letter codes found in polymer-specific config file.

    - **Recommended:** ["b", "y"] for peptides fragmented via CID.

- "ms2_adducts" :

    - **Description:** list of adducts to add to MS2 fragments other than intrinsic adducts specified in polymer-specific config file.

    - **Options:** any list of adducts found in GlobalChemicalConstants.

    - **Recommended:** ["Na"] for standard experiments in positive mode.

- "ms2_losses" :

    - **Description:** specifies whether to include side chain-specific neutral loss products for MS2 fragment ions.

    - **Options:** true or false.

    - **Recommended:** true.

- "ms2_max_neutral_losses" :

    - **Description:** specify upper cap on number of side chain-specific neutral loss products to include for MS2 fragments.

    - **Options:** null or integer.

    - **Recommended:** null.

- "ms2_loss_products_adducts" :

    - **Description:** specify whether to add any adducts on to MS2 neutral loss products.

    - **Options:** null or list of any adducts found in GlobalChemicalConstants.py.

    - **Recommended:** null.

- "add_signatures" :

    - **Description:** specify whether to add monomer-specific signature ions for MS2 fragments.

    - **Options:** true or false.

    - **Recommended:** true.

- "signatures" :

    - **Description:** if monomer MS2 signatures are to be added, this specifies any particular subsets of signatures that are to be added. If "add_signatures"==true and null is supplied, ALL signature types will be added.

    - **Options:** null or list of signature types specified in polymer-specific config file.

    - **Recommended:** null.

- "min_z" :

    - **Description:** specifies minimum charge of MS2 ions.

    - **Options:** integer.

    - **Recommended:** 1.

- "max_z" :

    - **Description:** specifies maximum charge of MS2 ions.

    - **Options:** integer.

    - **Recommended:** 1 for standard CID of small to medium-sized polymers, with the exception of polymers with excess of intrinsically charged monomers and / or large multi-metal centre transition metal complexes.

Extractor Parameters
====================

general extractor parameters
----------------------------

- "error" :

    - **Description:** specifies error tolerance threshold for matching target ions to ions observed in mass spectra. This can be supplied as an absolute value (amu) or relative error (i.e. ppm).

    - **Options:** float.

    - **Recommended:** 0.01 (absolute).

- "err_abs" :

    - **Description:** specifies whether error units are in absolute mass units or ppm.

    - **Options:** true or false.

    - **Recommended:** true.

- "min_ms2_peak_abundance" :

    - **Description:** specifies minimum relative abundance of most intense matching MS2 peak for a target sequence in observed MS2 spectra, expressed as a % of the most intense peak observed in the spectrum.

    - **Options:** float (range = 0 to 100).

    - **Recommended:** 90-100.

- "pre_run_filter" :

    - **Description:** specifies whether to filter spectra before beginning sequence screening. If true, spectra will be filtered through pre_screen_filters and sequences will be screened against truncated data set of spectra that pass these filters.

    - **Options:** true or false.

    - **Recommended:** true.


pre-screen filters
------------------

- "min_rt" :

    - **Description:** specifies minumum retention time of raw spectra. Any spectra with retention time lower than this will be discarded from consideration.

    - **Options:** float.

    - **Recommended:** 0 unless you know your chromatography well enough.

- "max_rt" :

    - **Description:** specifies maximum retention time of raw spectra. Any spectra with retention time higher than this will be discarded from consideration.

    - **Options:** float.

    - **Recommended:** run time (in minutes) of your method unless you know exactly what you are doing with this.

- "essential_signatures" :

    - **Description:** list of monomers that MUST have signatures present in MS2 spectra. Any MS2 spectra without these signatures will be discarded.

    - **Options:** null or list of monomer one-letter codes. These monomers MUST have signature ions in polymer-specific config file.

    - **Recommended:** null, unless you know for certain every sequence has these monomers AND  that these monomer signatures will always show up in MS2 spectra.

- "signature_types" :

    - **Description:** specifies list of signature codes to include in essential signatures. These MUST correspond to signature types in polymer-specific config file.

    - **Options:** any list of signature type found in polymer-specific config file.

    - **Recommended:** ["Im"] for peptides using depsipeptide config file, as this signature class represents immonium ions.

- "signature_ms_level" :

    - **Description:** specifies which MS level signatures will be present.

    - Options integer

    - **Recommended:** 2, as presently we only have data for MS2 signatures.

- "massdiff_bins" :

    - **Description:** *to be completed*

    - **Options:** true or false.

    - **Recommended:** *to be completed*

- "ms2_precursors" :

    - **Description:** list of precursor ions for matching to MS2 parents in observed MS2 spectra. Any MS2 spectra that do not have a parent matching one or more of these precursors will be discarded from consideration.

    - **Options:** null or list of floats.

    - **Recommended:** null unless doing very targeted screening for specific small subset of products (which general screen should pick up anyway...).

- "min_MS1_total_intensity" :

    - **Description:** specifies minimum total intensity of MS1 spectra. Any MS1 spectra with total intensity lower than this will be discarded from consideration.

    - **Options:** null or float.

    - **Recommended:** null.

- "min_MS2_total_intensity" :

    - **Description:** specifies minimum total intensity of MS2 spectra. Any MS2 spectra with total intensity lower than this will be discarded from consideration.

    - **Options:** null or float.

    - **Recommended:** null.

- "min_MS1_max_intensity" :

    - **Description:** specifies minimum maximum intensity of MS1 spectra. Any MS1 spectra with maximum intensity lower than this will be discarded from consideration.

    - **Options:** null or float.

    - **Recommended:** null.

- "min_MS2_max_intensity" :

    - **Description:** specifies minimum maximum intensity of MS2 spectra. Any MS2 spectra with maximum intensity lower than this will be discarded from consideration.

    - **Options:** null or float.

    - **Recommended:** null.

Postprocess Parameters
======================

- "exclude_frags" :

    - **Description:** list of specific fragment ids to exclude from consideration in confidence assignments ONLY (will still be included in retention time assignment).

    - **Options:** null or any list of valid fragment ids.

    - **Recommended:** null.

- "optional_core_frags" :

    - **Description:** list of specific fragment ids to exclude from consideration in confidence assignments IF they are absent.

    - Options null or any list of valid fragment ids.

    - **Recommended:** ["b1"] for peptides in positive mode, as this particular MS2 fragment is notorious for not showing up in MS2 spectra of peptides.

- "core_linear_series" :

    - **Description:** list of fragment series one letter codes for fragment series that should be used in confidence calculations. NOTE: it is assumed that each series in this list will receive an equal weighting in confidence calculation; if this is not suitable for your chemistry, please let us know!!!

    - **Options:** list of valid fragment one letter codes. NOTE: this MUST not include signatures.

    - **Recommended:** ["b", "y"] for peptides in positive mode (["y"] for peptides in negative mode).

- "excluded_fragments" :

    - **Description:** list of specific fragment ids to exclude from confidence assignments AND retention time assignments.

    - **Options:** null or any list of specific fragment ids.

    - **Recommended:** null.

- "dominant_signature_cap" :

    - **Description:** specifies an upper cap on confidence assignments for sequences that are expected to have abundant monomer-specific MS2 signature(s) but have one or more of these missing.

    - **Options:** float (range = 0 to 100).

    - **Recommended:** 70.

- "essential_fragments" :

    - **Description:** specifies specific fragment ids that MUST be present for a sequence to be assigned with confidence.

    - **Options:** null or list of specific valid fragment ids.

    - **Recommended:** null.

- "subsequence_weight" :

    - **Description:** specifies weighting to place on continuous fragment coverage in calculating overall confidence of a sequence assignment. This is supplied as a decimal fraction.

    - **Options:** float (range = 0 to 100).

    - **Recommended:** we are still trying to work this out. Probably between 0.25 and 0.5.

- "min_rt" :

    - **Description:** mininimum retention time to assign to confirmed sequences.

    - **Options:** null or float.

    - **Recommended:** null.

- "max_rt" :

    - **Description:** maximum retention time to assign confirmed sequences.

    - **Options:** null or float.

    - **Recommended:** null.

- "Rt_bin" :

    - **Description:** specifies minimum separation of MS1 peaks in extracted ion chromatograms (units = minutes).

    - **Options:** float.

    - **Recommended:** 0.25 for LC-MS, 0 without chromatography.

- "backup_Rt_bin" :

    - **Description:** specifies a back-up retention time bin if peaks cannot be assigned with Rt_bin. This may be useful in cases where particular sequence subsets do not separate well using chosen chromatography method (which happens sometimes for complex mixtures).

    - **Options:** float (should be lower than Rt_bin).

    - **Recommended:** 0.1.

- "ms2_Rt_bin" :

    - **Description:** specifies maximum discrepancy between retention time of unique MS2 fragment and MS1 precursor.

    - **Options:** float.

    - **Recommended:** 0.5 (this may be too generous).

- "ms2_Rt_flexible" :

    - **Description:** specifies whether to allow for the possibility of a mismatch between observed MS2 retention time and MS1 precursor within specified ms2_Rt_bin. If true, ms2_Rt_bin tolerance will be widened until a retention time can be assigned.

    - **Options:** true or false.

    - **Recommended:** we need to test this urgently.

- "min_viable_confidence" :

    - **Description:** specifies minimum confidence score for sequences to be assigned final retention time and intensity values. Sequences with scores lower than this will not be processed for retention time and / or intensities.

    - **Options:** float.

    - **Recommended:** 55-60.

- "min_relative_intensity" :

    - **Description:** specifies minimum relative intensity (as % intensity of most intense peak) in EICs. Any peaks with relative intensities lower than this will be discarded, which may lead to less abundant sequences not being assigned.

    - **Options:** null or float (range = 0 to 100).

    - **Recommended:** null, unless you are deliberately attempting to assign retention time of very low confidence sequences (which probably won't work anyway..).

- "plot EICs" :

    - **Description:** option to plot all EICs and save the plots as png files into your output folder. NOTE: not finished yet.

    - **Options:** true or false.

Directories
===========

- "ripper_folder" :

    - **Description:** specifies location of folder containing mzml ripper JSON files.

    - **Options:** string file path to folder.

    - **Recommended:**

- "output_folder" :

    - **Description:** specifies location of output folder, where output data will be saved.

    - **Options:** full string file path to output folder.

    - **Recommended:** scapa data folder (obviously).