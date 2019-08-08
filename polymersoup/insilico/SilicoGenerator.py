"""
This file contains functions for passing on silico parameters to insilico
functions for generating theoretical MS and MSMS sequence data

"""

from .MS1_silico import *
from .MS2_silico import *
from ..parameter_handlers import * 
from .silico_helpers.insilico_helpers import * 



def generate_MSMS_sequence_dict(
    parameters_file="C:/Users/group/PolymerMassSpec/Examples/InputParams_Test.json",
    write=True,
    output_folder=None,
    ms1=True,
    ms2=True
):
    """
    This function is used to generate a full MS and MSMS in silico sequence
    dict in format: {
                        seq: {
                            'MS1' : [m/z...],
                            'MS2': {
                                'frag': [m/z...],
                                'signatures': {
                                    'signature': [m/z...]
                                },
                                'unique_fragments' ['unique_frag'...]
                            }
                        }
        where seq = sequence string, m/z = m/z value of a given ion, 'frag' =
        MS2 fragment id string, 'signature' = MS2 fragment id string for a
        monomer-specific signature fragment, 'unique_frag' = MS2 fragment id
        string for a fragment with m/z unique to that sequence (i.e. not found
        in MS2 fragments of other, isobaric sequences)

    Args:
        parameters_file (str, optional): full file path to input parameters
            .json file. Defaults to "C:/Users/group/PolymerMassSpec/Examples/InputParams_Test.json".
        write (bool, optional): specifies whether to write sequence dict to
            file. Defaults to True.
        output_folder (str, optional): full file path to output folder to put
            output file; if None, file is dumped in same directory as
            parameters_file. Defaults to None.
        ms1 (bool, optional): specifies whether to generate theoretical MS1
            data for sequences. Defaults to True.
        ms2 (bool, optional): specifies whether to generate theoretical MS2
            data for sequences. Defaults to True.

    Returns:
        full_MSMS_dict (dict): MSMS sequence dict, the format of which
            is described in detail above
    """
    # reads in silico parameters from input parameters json
    silico_params = return_silico_parameters(parameters_file)

    if ms1:
        start = time.time()
        ms1 = silico_params["MS1"]

        MS1 = generate_ms1_mass_dictionary_adducts_losses(
            ms1["monomers"],
            ms1["max_length"],
            ms1["ms1_adducts"],
            silico_params["mode"],
            ms1["min_z"],
            ms1["max_z"],
            ms1["losses"],
            ms1["max_neutral_losses"],
            ms1["loss_products_adducts"],
            ms1["min_length"],
            ms1["chain_terminators"],
            ms1["universal_rxn"],
            ms1["terminal_tags"]["0"],
            ms1["terminal_tags"]["-1"]
        )
        end = time.time()
        time_elapsed = end-start
        print(f'for {len(MS1)} sequences, time taken for MS1 = {time_elapsed} seconds')
        print(f'time taken per sequence for MS1 = {time_elapsed/len(MS1)} seconds')

    start = time.time()
    ms2 = silico_params["MS2"]
    MS2 = generate_ms2_mass_dictionary(
            MS1.keys(),
            ms2["fragment_series"],
            ms2["ms2_adducts"],
            silico_params["mode"],
            ms2["add_signatures"],
            ms2["signatures"],
            ms2["ms2_losses"],
            ms2["ms2_max_neutral_losses"],
            ms2["ms2_loss_products_adducts"],
            ms2["uniques"]
    )
    end = time.time()
    time_elapsed = end-start
    print(f'for {len(MS2)} sequences, time taken for MS2 = {time_elapsed} seconds')
    print(f'time taken per sequence for MS2 = {time_elapsed/len(MS2)} seconds')
    full_MSMS_dict = {}

    for sequence in MS1:
        full_MSMS_dict[sequence] = {
            'MS1': MS1[sequence],
            'MS2': MS2[sequence]
            }

    full_MSMS_dict = add_peak_lists_massdict(full_MSMS_dict)

    if write:
        if not output_folder:
            output_folder = os.path.dirname(parameters_file)

        write_file = generate_insilico_writefile_string(
            output_folder,
            silico_params
        )

        write_to_json(full_MSMS_dict, write_file)

    return full_MSMS_dict

generate_MSMS_sequence_dict()

def generate_MS1_compositional_dict(
    silico_parameters
):

    # check if silico_parameters is provided as a file path to input parameters
    # .json file or already opened as a dict; if file path, open as dict
    if type(silico_parameters) == str:
        silico_parameters = open_json(silico_parameters)['silico_parameters']

    # load starting monomers
    monomers = silico_parameters['MS1']['monomers']

    # load minimum and maximum sequence length
    min_length = silico_parameters['MS1']['min_length']
    max_length = silico_parameters['MS1']['max_length']

    # load mass spec mode - either pos or neg
    mode = silico_parameters['mode']

    # load MS1 adducts
    adducts = silico_parameters['MS1']['ms1_adducts']

    # load minimum and maximum MS1 charge state for sequences
    min_z = silico_parameters['MS1']['min_z']
    max_z = silico_parameters['MS1']['max_z']

    # check whether side chain-specific neutral loss products are to be
    # included in MS1 compositional info, the maximum number of such losses
    # allowed per sequence, and any associated adducts
    losses = silico_parameters['MS1']['losses']
    max_neutral_losses = silico_parameters['MS1']['max_neutral_losses']
    loss_products_adducts = silico_parameters['MS1']['loss_products_adducts']

    # find out whether all monomers are universally cross-reactive
    universal_rxn = silico_parameters['MS1']['universal_rxn']

    # list any monomers that terminate elongation
    chain_terminators = silico_parameters['MS1']['chain_terminators']

    # retrieve any tags to add on to termini
    terminal_tags = silico_parameters['MS1']['terminal_tags']
    start_tags = terminal_tags['0']
    end_tags = terminal_tags['-1']

    MS1_dict = generate_ms1_mass_dictionary_adducts_losses(
        monomers,
        max_length,
        adducts,
        mode,
        min_z,
        max_z,
        losses,
        max_neutral_losses,
        loss_products_adducts,
        min_length,
        chain_terminators,
        start_tags,
        end_tags,
        sequencing=False
    )

    return MS1_dict
