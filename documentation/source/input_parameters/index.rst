.. oligoss documentation master file, created by
   sphinx-quickstart on Tue Sep  8 17:24:36 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. _Input-Parameters:

#####################
Input Parameters File
#####################

Input parameter / experiment run files define all parameters that may change between
individual sequencing experiments. Many of these parameters have default values, some of
which are determined by instrument or oligomer class. The input parameter file is in JSON
format, and is divided into four sections:

* :ref:`core parameters<Core-Parameters>`
* :ref:`silico parameters<Silico-Parameters>`: parameters relevant to *in silico* fragmentation and generation of theoretical precursor-MS2 product ion libraries.
* :ref:`extractor parameters<Extractor-Parameters>`: parameters relevant to screening, filtering and extracting data from observed MS and MS2 spectra.
* :ref:`postprocess parameters<Postprocess-Parameters>`: parameters relevant to assigning final sequence assignments and confidence scores to extracted data.

.. toctree::
   :maxdepth: 4

   general_parameters
   silico_parameters
   extractor_parameters
   postprocess_parameters

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
