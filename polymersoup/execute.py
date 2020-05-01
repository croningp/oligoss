"""
PolymerSoup version 2.0

This file is used to run all PolymerSoup experiments.
Artificial Life Team, Cronin Group, 2020
"""
import os
import sys
import json
from .utils.parameter_handlers.parameter_handlers import generate_parameters

def run_polymersoup(input_file):
    """
    Reads input file and runs a PolymerSoup workflow

    Args:
        input_file (str): path to input file
    """

    #  open input parameters file as dict
    params_json = os.path.abspath(input_file)
    with open(params_json) as f:
        params_dict = json.load(f)

    #  generate instance of Parameters class from input parameters dict
    run_params = generate_parameters(
        params_dict=params_dict)

    #  temp print of run parameters for debugging
    print(vars(run_params))


run_polymersoup(input_file=sys.argv[1])
