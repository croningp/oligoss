.. _Installation-Requirements:

#############################
Installation and Requirements
#############################

.. _Dependencies:

Dependencies
============

**Python (version 3.6.0 or later)**

OLIGOSS was written in Python 3.7.5, but should be compatible with version 3.6
or later.

**mzML**

* Mass spectrometry data must be in .mzML file format. mzML files can be generated from a variety of vendors, including Proteowizard MS Convert, which is freely available here: http://proteowizard.sourceforge.net/download.html

* mzMLs are automatically converted to JSON format using Graham Keenan's mzmlripper, documentation for which can be found on `GitHub <https://github.com/croningp/mzmlripper.git>`_.

.. _System-Requirements:

System Requirments
==================

OLIGOSS was developed and tested on Ubuntu 19.10, and should therefore be
compatible with any Unix OS. As of version 0.0.3, OLIGOSS is incompabitible
with Windows. Windows compatibility will be introduced in a later version.
Luckily for Windows 10 users, the latest distributions of Ubuntu can be easily
installed and run from Windows 10. For instructions on how to set up Ubuntu on Windows 10
and install OLIGOSS, see our instructions `here<Install-Windows>`.

.. _Standard-Install:

Standard Installation
=====================

OLIGOSS is avaialable through The Python Package Index (`PyPI <https://pypi.org/project/oligoss/>`_).

.. _Install-PyPi:

Install From PyPI
-----------------

To install OLIGOSS from PyPi, use pip::

   pip3 install oligoss

.. _Install-Source:

Download From GitHub
--------------------

Alternatively, the source code can be cloned directly from the `OLIGOSS GitHub repository <https://github.com/croningp/oligoss.git>`_.

.. _Install-Windows:

Windows 10 Installation
-----------------------

OLIGOSS was largely written and tested using Python 3.7.5 on operating system Ubuntu 19.04. It is currently not compatible with
Windows. Luckily for Windows 10 users, Ubuntu can be easily installed and run for Windows computers using the Windows Subsystem for Linux 
(WSL).

Installing WSL and Ubuntu on Windows 10
***************************************

#. To install WSL and Ubuntu for Windows 10, follow the instructions from `Microsoft <https://docs.microsoft.com/en-us/windows/wsl/install-win10>`
   * For new users, we recommend using the recommended options in the Microsoft instructions (WSL2, most recent distribution of Ubuntu).

#. Once Ubuntu has been installed succesfully, launch Ubuntu.

#. Within the Ubuntu shell, make sure the package list is up-to-date with the following command::
   
   sudo apt update

#. Install `pip <https://pip.pypa.io/en/stable/installing/>` using the command::
   
   sudo apt install python3-pip

#. Finally, use pip to install OLIGOSS::
   
   pip3 install oligoss

.. _Running-OLIGOSS:

Running OLIGOSS
===============

To run an OLIGOSS sequencing workflow, run the following command::

    python3 -m oligoss -i input_params.json -r ripper_folder -o output_folder

**input_params.json = input parameters file**

    This should contain all relevant input parameters for executing an
    OLIGOSS sequencing workflow (see :ref:`Input-Parameters`).

**ripper_folder = data directory**.

    This argument can either be passed in via the command line directly (as above) or specified in the input parameters file using the data_folder parameter.
    This folder should contain input MS data in either mzML or ripper JSON format.

**out_dir = output directory**.

    This argument can either be passed in via the command line directly (as above) or specified in the input parameters file using the output_folder parameter.
    All output data will be dumped to this folder.

