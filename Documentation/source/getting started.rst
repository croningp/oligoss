Getting Started
###############

Dependencies
============

**Python (version 3.6.0 or later)**

- This package was written in Python 3.7, but should be compatible with version 3.6 or later.

**mzML**

- Mass spectrometry data must be in .mzML file format. mzML files can be generated from a variety of vendors, including Proteowizard MS Convert, which is freely available here: http://proteowizard.sourceforge.net/download.html

**mzmlripper**

- mzML files must be converted to JSON format using Graham Keenan's mzmlripper, instructions and documentation for which can be found here: http://datalore.chem.gla.ac.uk/Origins/mzmlripper.git

Installation
============

The polymermassspec repo can be cloned from:
http://datalore.chem.gla.ac.uk/DBD/polymermassspec.git

Running PolymerMassSpec
=======================

(1) Convert raw data files to JSON files using mzmlripper.

(2) Ensure Input Parameters JSON file is filled in with correct run parameters (see "Run Parameters" section).

(3) Run executable script:

    - Experiments are run from the executable, which can be found here: polymersoup/executable.py
    - To run the script: **python -m polymersoup.executable input_parameters_file.json**
