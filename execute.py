"""
PolymerSoup version 2.0

This file is used to run all PolymerSoup experiments.
Artificial Life Team, Cronin Group, 2020
"""
import os
import time
import logging

from logging.config import dictConfig

from polymersoup.utils.logging_utils.logger_utils import (
    generate_logger_config_dict,
    get_git_info_for_log
)
from polymersoup.utils.parameter_handlers.parameter_handlers import (
    generate_parameters
)
from polymersoup.utils.run_utils import arg_parser
from polymersoup.workflows.workflows import (
    exhaustive_screen_multi
)

#  get Namespace object from CLI
args = arg_parser()

def run_polymersoup(input_file, data_folder=None, out_folder=None):
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

    if not data_folder:
        data_folder = run_params.data_folder
        if not data_folder:
            raise Exception(run_params.data_folder)
    if not out_folder:
        out_folder = run_params.output_folder
        if not out_folder:
            raise Exception(
                "output directory must be specified in input parameters or\n"
                "CLI. For help, run python execute.py -h")
        if not os.path.exists(out_folder):
            os.mkdir(out_folder)

    #  set log file for root logger
    log_file = os.path.join(out_folder, f"run_log_{time.time()}.log")
    if not os.path.exists(os.path.dirname(log_file)):
        os.mkdir(os.path.dirname(log_file))
    dictConfig(generate_logger_config_dict(log_file))

    #  log git user, version and branch
    logging.info(get_git_info_for_log())
    #  log message at start of screen
    logging.info("beginning exhaustive screen")

    exhaustive_screen_multi(
        params=run_params,
        data_folder=data_folder,
        out_folder=out_folder
    )


if __name__ == '__main__':
    run_polymersoup(
        input_file=args.input,
        data_folder=args.ripper,
        out_folder=args.output
    )
