import sys
import logging
import argparse
import multiprocessing

def set_max_cores(params: object) -> int:
    """
    Sets maximum cores to use for multiprocessing based on number of cores to
    be left available.

    Args:
        params (Parameters): parameters object.

    Returns:
        int: maximum number of cores to use.
    """

    #  get max cores available
    max_cores = multiprocessing.cpu_count()

    #  check whether upper limit has been specified in input params
    if params.max_cores:
        max_cores = min(max_cores, params.max_cores)

    #  if cores need to be left avaialble for other use, further reduce
    #  maximum used in workflows
    if params.free_cores:
        max_cores = max_cores - params.free_cores

    #  log max_cores
    logging.info(
        f"maximum numnber of cores available for multiprocessing = {max_cores}")

    #  set all available cores for use
    return max_cores

def check_child_process_fork() -> bool:
    """
    Checks OS to determine whether child process default start method is
    fork or spawn.

    Returns:
        bool: True or False
    """
    log_msg = "OS = "
    if sys.platform.startswith("linux"):
        logging.info(f"{log_msg}{sys.platform}. Child processes forked.")
        return True
    logging.info(f"{log_msg}{sys.platform}. Child processes spawned")
    return False

def arg_parser() -> object:
    """
    Sets arguments for Polymersoup execute.py CLI.

    Returns:
        object: Namespace object with args as properties.
    """

    #  create ArgumentParser object
    my_parser = argparse.ArgumentParser(
        description="executes a Polymersoup sequencing workflow",
        epilog="Happy Sequencing!"
    )

    #  add arguments for input (run parameters file), ripper directory (folder
    #  containing ripper JSONs and / or mzML files) and output directory for
    #  dumping results
    my_parser.add_argument(
        "-input",
        action="store",
        required=True,
        help="full filepath to input run parameters file"
    )
    my_parser.add_argument(
        "-ripper",
        action="store",
        required=False,
        help="directory containing mzML files or ripper JSONs. If this is not\n"
        "specified in command line, it must be supplied in input\n"
        "parameters"
    )
    my_parser.add_argument(
        "-output",
        action="store",
        required=False,
        help="directory containing mzML files or ripper JSONs. If this is not\n"
        "specified in command line, it must be supplied in input parameters"
    )

    return my_parser.parse_args()
