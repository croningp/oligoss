.. _Core-Parameters:

Core Parameters
===============

Input parameter files must be in JSON format.
The following core parameters must be explictly defined in each input parameters file:

#. **mode**:
    *  Description: specifies the charge of ions.
    *  Type: `str`
    *  Options: either "pos" or "neg" for positive and negative ion mode, respectively.


#. **monomers**:
    * Description: monomers used in experiment.
    * Type: `List[str]`
    * Options: any combination of monomer one letter codes from polymer-specific configuration files.


#. **screening method**:
    * Description: specifies which sequencing workflow to use.
    * Type: `str`
    * Options: currently only option is "exhaustive".


#. **polymer_class**:
    * Description: specifies the backbone of oligomers being fragmented.
    * Type: `str`
    * Options: any valid alias associated with a polymer config file (see **Polymer Configs**).

In addition to these required parameters, optional instrument configuration files can be defined.
**NOTE**: if these are defined, some parameters for silico, extractors and postprocessing are
assigned default instrument values.

#. **instrument**:
    * Description: specifies the mass spec model used to acquire MS/MS data.
    * Type: `str`, optional
    * Options: any valid alias for a pre-configured instrument file.


#. **chromatography**:
    * Description: specifies the chromatography used to separate products prior to detection via MS.
    * Type: `str`, optional
    * Options: any valid alias for a pre-configured chromatography file.

#. **max_cores**:
   * Description: specifies the maximum number of logical cores to utilise for sequecing workflows -
      OLIGOSS utilises multiprocessing to speed up data processing.
   * Type: `int`, optional
   * Options: any number of cores. See Python's `documentation for multiprocessing <https://docs.python.org/3.7/library/multiprocessing.html>`_ for more information.
   * Default: 4.

#. **free_cores**:
     * Description: specifies the minimum number of free logical cores to be left unused during sequencing workflows.
     * Type: `int`, optional
     * Options any number of cores from 0 to N - 1, where N = number of logical cores (threads) in system. 
     * Default: 2.

#. **hard_memory_limit**:
    * Description: maximum amount of virtual memory to allocate to the OLIGOSS run (in Gb).
    * Type: `float`, optional
    * Options: RAM usage (in Gb). If a run exceeds this limit, it will be terminated prematurely.
    * Default: None.

#. **relative_memory_limit**:
   * Description: maximum amount of virtual memory to allocate to OLIGOSS run (as % of total virtual memory).
   * Type: `float`, optional
   * Options: any value > 0 and <= 100.
   * Default: 90.