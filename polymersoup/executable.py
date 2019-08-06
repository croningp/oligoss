"""
This file will be populated with high level functions used to perform
full sequencing runs, calling upon functions from the following modules:

1. extractors
2. insilico
3. postprocessing
4. filehandler
"""
from .extractors.sequence_screening import *
from .insilico.SilicoGenerator import *
from .filehandler import *

def launch_screen(input_parameters_file):
    """
    This function reads an input parameters .json file and decides what
    kind of screening method to use, then performs sequencing screen

    Args:
        input_parameters_file (str): full file path to input parameters .json
            file
    """
    # load parameters
    parameters_dict = open_json(input_parameters_file)

    if parameters_dict['screening_method'] == 'exhaustive':
        exhaustive_screen(parameters_dict)


def exhaustive_screen(
    parameters_dict
):
    """
    This function performs an exhaustive or 'brute force' screen for all
    possible sequences arising from input monomers and constraints set, read
    directly from input parameters .json file

    Args:
        parameters_dict (dict): input parameters dictionary in standard
            input parameters format
    """

    # load parameters for in silico operations
    silico_params = parameters_dict['silico_parameters']

    # load parameters for data extraction operations
    extractor_params = parameters_dict['extractor_parameters']

    directories = parameters_dict['directories']

    ripper_folder = directories['ripper_folder']

    # generate compositional silico dict for screening MS1 EICs
    compositional_silico_dict = generate_MS1_compositional_dict(silico_params)

    # filter mass spectra to remove useless data before screening
    pre_filter_rippers(
        ripper_folder,
        parameters_dict
    )

    # retrieve file paths to fitlered mass spectra
    filtered_ripper_folder = os.path.join(ripper_folder, 'filtered')
    filtered_rippers = [
        os.path.join(filtered_ripper_folder, file)
        for file in os.listdir(filtered_ripper_folder)
    ]

    MS1_EICs = extract_MS1
