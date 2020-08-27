import sys
import time
import logging
import argparse
import multiprocessing

class RunInfo():
    """
    Class for keeping track of important system info - operating system, number
    of threads available for use.
    """

    def __init__(self, params):
        self._params = params
        self._total_cores = multiprocessing.cpu_count()
        self.max_cores = self.set_max_cores()
        self._os = sys.platform
        if self._os.startswith("linux"):
            self.init_method = "forked"
        else:
            self.init_method = "spawned"

    def set_max_cores(self):
        max_cores = self._total_cores
        if self._params.max_cores:
            max_cores = min(self._params.max_cores, max_cores)
        if self._params.free_cores:
            if max_cores > self._total_cores - self._params.max_cores:
                max_cores = self._total_cores - self._params.max_cores
        if max_cores < 1:
            raise Exception(f"max cores must be > 0, not {self.max_cores}")
        return max_cores

    def __str__(self):
        """
        str method to be used for logging and debugging.

        Returns:
            str: message with info on number of cores available for workflows
                and the operating system.
        """

        #  generate message for number of cores avaialble.
        core_log = (
            f"Total cores on computer = {self._total_cores}. Max available\n"
            f"for use = {self.max_cores}. Number of free cores specified\n"
            f"in input parameters = {self._params.free_cores}.\n")
        if self.max_cores == 1:
            core_log += (
                "Only 1 core is avaiable. Polymersoup performs best when it\n"
                "can use multiple cores to execute code. Please check input\n"
                "parameters for maximum allocated cores and number of free /\n"
                "reserved cores.\n")

        #  generate message for OS and how this affects the workflow
        os_log = f"Operating system = {self._os}"
        if self._os.startswith("linux"):
            os_log += (
                "forked exhaustive screen will be run without using MongoDB")
        else:
            os_log += (
                "spawned exhaustive screen will be run with MongoDB")

        return core_log + os_log

    def log_record(self):
        logging.info(str(self))

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

def log_screening_run(screening_func):
    """
    Very basic logging decorator for screening functions that are used to screen
    a single ripper.

    Args:
        screening_func (function): screening function with args and / or kwargs.
    """

    def log_screen(*args, **kwargs):
        """
        Attempts screening function, logs total run time and any runtime errors.

        Returns:
            Optional[Any]: either NoneType in the event of failure or whatever
                value screening_func returns.
        """

        #  keep track of which ripper file is being screened
        if kwargs:
            if "ripper_name" in kwargs:
                ripper_file = kwargs["ripper_name"]
            elif "ripper" in kwargs:
                ripper_file = kwargs["ripper"]
        else:
            if args:
                ripper_file = args[0]
            else:
                ripper_file = "unknown"

        #  retrieve basic run information - i.e. threads used, init method for
        #  child processes
        if kwargs and "run_info" in kwargs:
            run_info = kwargs["run_info"]
        else:
            run_info = [arg for arg in args if type(arg) == RunInfo]
            if run_info:
                run_info = run_info[0]

        #  assume result is NoneType. This will be updated to screening_func
        #  return value if screen is succesful
        result = None

        #  attempt the screen, keeping track of time taken and log both success
        #  and failure
        try:
            start = time.time()
            logging.info(
                f"attempting to screen {ripper_file} using"
                f"{screening_func.__name__}")
            end = (time.time() - start) / 60
            result = screening_func(*args, **kwargs)
            logging.info(
                f"screening of {ripper_file} successful after {end} minutes.")
        except Exception as e:
            end = (time.time() - start) / 60
            logging.exception(
                f"screening of {ripper_file} failed after {end} minutes "
                f"with exception: {e}")
        finally:
            if run_info:
                logging.info(
                    f"screen run on {run_info.max_cores} threads with "
                    f"{run_info.init_method} child processes")
        return result

    return log_screen
