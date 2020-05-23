"""
PolymerSoup version 2.0

This file is used to run all PolymerSoup experiments.
Artificial Life Team, Cronin Group, 2020
"""
import os
import sys
from polymersoup.utils.parameter_handlers.parameter_handlers import (
    generate_parameters
)
from polymersoup.workflows.workflows import exhaustive_screen

def run_polymersoup(input_file, data_folder, out_folder):
    """
    Reads input file and runs a PolymerSoup workflow

    Args:
        input_file (str): path to input file
    """

    #  open input parameters file as dict
    params_json = os.path.abspath(input_file)

    #  generate instance of Parameters class from input parameters dict
    run_params = generate_parameters(
        params_json=params_json)

    exhaustive_screen(
        params=run_params,
        data_folder=data_folder,
        out_folder=out_folder
    )


run_polymersoup(
    input_file=sys.argv[1],
    data_folder=sys.argv[2],
    out_folder=sys.argv[3]
)
