"""
This file contains functions for generating theoretical MS and MSMS data for
ALife polymer chemistry.
"""
from .MS1_silico import *
from .MS2_silico import *
from .. parameter_handlers import *
def __init__():
    """
    Initialises SilicoGenerator and imports constants from polymer config file

    """

def generate_MSMS_sequence_dict(
    parameters_file="C:/Users/group/PolymerMassSpec/Examples/InputParams_Test.json",
    write=True,
    output_folder=None,
    ms1=True,
    ms2=True
):

    # reads in silico parameters from input parameters json
    silico_params = return_silico_parameters(parameters_file)

    if ms1:
        start = time.time()
        ms1 = params_dict["MS1"]

        MS1 = generate_ms1_mass_dictionary_adducts_losses(
            ms1["monomers"],
            ms1["max_length"],
            ms1["ms1_adducts"],
            ms1["mode"],
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
    print(f'doing in silico fragmentation for {len(MS1)} sequences')
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
