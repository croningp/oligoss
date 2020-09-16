Getting Started
###############

Dependencies
============

**Python (version 3.6.0 or later)**

This package was written in Python 3.7, but should be compatible with version 3.6 or later.

**mzML**

* Mass spectrometry data must be in .mzML file format. mzML files can be generated from a variety of vendors, including Proteowizard MS Convert, which is freely available here: http://proteowizard.sourceforge.net/download.html
* mzML's are automatically converted to JSON format using Graham Keenan's mzmlripper, documentation for which can be found here: http://datalore.chem.gla.ac.uk/Origins/mzmlripper.git

Installation
============

The polymermassspec repo can be cloned from:
http://datalore.chem.gla.ac.uk/DBD/polymermassspec.git

Packages Required
=================

A list of packages needed to run polymersoup can be found in requirements.txt

Running the following command within the directory, after cloning the repository, will install all necessary packages::

    pip install -r requirements.txt 


Running PolymerMassSpec
=======================

To run a Polymersoup sequencing workflow, run the following command::

    python polymersoup/execute.py -i input_params.json -r ripper_folder -o output_folder

**input_params.json = input parameters file**

    This should contain all relevant input parameters for executing a Polymersoup sequencing workflow (see Input Parameters, below).

**ripper_folder = data directory**.

    This argument can either be passed in via the command line directly (as above) or specified in the input parameters file using the data_folder parameter.
    This folder should contain input MS data in either mzML or ripper JSON format.

**out_dir = output directory**.

    This argument can either be passed in via the command line directly (as above) or specified in the input parameters file using the output_folder parameter.
    All output data will be dumped to this folder.

Sequencing Workflows
====================

There is currently only one sequencing workflow available in Polymersoup.

The **exhaustive screening** workflow is to be used for sequencing oligomers with well-characterised fragmentation pathways and known monomer libraries.

For oligomer classes with poorly characterised fragmentation pathways and / or data with unknown monomers, knew workflows based on previous **mass difference** screen and Kendrick Analysis will be coming shortly.

.. image:: img/exhaustive_workflow_snip.png
    :width: 600
    :align: center
