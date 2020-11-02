import os
import sys
import time
import signal
import socket
import psutil
import pickle
import logging
import argparse
import multiprocessing

from functools import wraps

if sys.platform.startswith("linux"):
    import resource

class RunInfo():
    """
    Contains important info on system used for run - cores, available memory
    etc..

    Args:
        params (Parameters): Parameters object.
    """

    def __init__(self, params):

        #  protected attributes:
        #  Parameters object, number of cores avaialable on machine, operating
        #  system and total RAM (excluding swap)
        self._params = params
        self._total_cores = multiprocessing.cpu_count()
        self._os = sys.platform
        self._total_memory = psutil.virtual_memory().total
        self._version = sys.version_info

        #  get init_method for child processes in multiprocessed workflows
        if self._os.startswith("linux"):
            self.init_method = "forked"
        else:
            self.init_method = "spawned"

        #  get local user account and host computer for run
        self.user = psutil.Process().username
        self.host = socket.gethostname()

        #  init warning message for logging any warnings for run
        self.warning_msg = ""

        #  get memory limits from input parameters and apply (if possible - this
        #  depends on operating system). If applied successfully, MemoryError
        #  will be thrown if total RAM usage of parent + child processes exceeds
        #  limit.
        self.memory_cap = self.set_memory_cap()
        self.apply_memory_cap()

        #  get process id for parent process. This will be used later if
        #  if required to terminate parent and child processes
        self.parent_id = os.getpid()

    @property
    def max_cores(self):
        """
        Sets maximum number of cores to use in run.

        Raises:
            Exception: raised if max_cores < 1

        Returns:
            int: maximum number of cores to use in multiprocessed runs.
        """

        max_cores = self._total_cores
        if self._params.max_cores:
            max_cores = min(self._params.max_cores, max_cores)
        if self._params.free_cores:
            max_cores = min(
                max_cores,
                self._total_cores - self._params.free_cores)
        if max_cores < 1:
            raise Exception(f"max cores must be > 0, not {max_cores}")
        return max_cores

    def __str__(self):
        """
        str method to be used for logging and debugging.

        Returns:
            str: message with info on number of cores available for workflows
                and the operating system.
        """
        #  generate message for host as well as user account used to run
        #  Oligomersoup
        user_log = (
            f"Running on host {self.host} with {self.user} as process owner.\n")

        version_log = f"Run using Python version: {self._version}.\n"

        #  generate message for number of cores avaialble.
        core_log = (
            f"Total cores on computer = {self._total_cores}. Max available"
            f"for use = {self.max_cores}. Number of free cores specified"
            f"in input parameters = {self._params.free_cores}.\n")
        if self.max_cores == 1:
            core_log += (
                "Only 1 core is avaiable. Oligomersoup performs best when it"
                "can use multiple cores to execute code. Please check input"
                "parameters for maximum allocated cores and number of free /"
                "reserved cores.\n")

        #  generate message for OS and how this affects the workflow
        os_log = f"Operating system = {self._os}"
        if self._os.startswith("linux"):
            os_log += (
                "forked exhaustive screen will be run without using MongoDB")
        else:
            os_log += (
                "spawned exhaustive screen will be run with MongoDB")

        #  generate message for available virtual memory and limits placed on
        #  memory usage from input parameters
        memory_log = (
            f"Total system memory = {self._total_memory} bytes"
            f"({self._total_memory/1E9} Gb). Memory limits specified from input"
            f"parameters = {self._params.relative_memory_limit * 100} %"
            "(relative)"
        )
        if self._params.hard_memory_limit:
            memory_log += (
                f"; and {self._params.hard_memory_limit} Gb (absolute).")
        else:
            memory_log += ". No absolute memory limit specified."

        #  final message includes info on user, cores, OS and virtual memory
        return f"{user_log}{core_log}{os_log}{memory_log}{version_log}"

    def log_record(self):
        logging.info(str(self))
        if self.warning_msg:
            logging.warning(self.warning_msg)

    def set_memory_cap(self):
        """
        Checks for memory caps (absolute memory cap in Gb or relative memory
        cap). Ensures sequencing run do not exceed memory allocation.

        Raises:
            Exception: raised if relative_memory_limit (passed in from input
                parameters is > 1)

        Returns:
            float: memory limit for sequencing workflow to use (in bytes)
        """
        #  check for relative memory limit - this is a decimal fraction of total
        #  available memory
        if (
            self._params.relative_memory_limit > 1
                or self._params.relative_memory_limit < 0):
            raise Exception(
                "relative memory limit must be expressed as decimal"
                "fraction (i.e. between 0 and 1)")
        mem_limit = self._total_memory * self._params.relative_memory_limit

        #  check for hard memory limit specified in parameters. This should be
        #  in units of Gb
        if self._params.hard_memory_limit:

            #  set memory limit as whichever is lower => relative or hard memory
            #  limit. Make sure hard memory limit is converted from Gb to bytes
            mem_limit = min(mem_limit, self._params.hard_memory_limit * 1E9)

        return mem_limit

    def apply_memory_cap(self):
        """
        Enforces limitations on total RAM usage for processes in Oligomersoup
        run. NOTE: currently only works for Linux.
        """

        #  depending on operating system, enforce an upper limit on RAM usage
        if self._os.startswith("linux"):

            #  set memory cap to avoid using too much RAM
            resource.setrlimit(
                resource.RLIMIT_AS, (self.memory_cap, self.memory_cap))

            #  set warning message to include info on enforced memory cap
            self.warning_msg += (
                "Total virtual memory usage for process and children capped at"
                f"{round(self.memory_cap / 1E9, 2)} Gb")
        else:
            #  set warning message to include info on lack of enforcable mem cap
            self.warning_msg += (
                "We are unable to enforce a memory limit as there is currently "
                f"no functionality for capping memory for {self._os} "
                "operating system.")

    @property
    def children(self):
        """
        Returns list of children processes.

        Returns:
            List[psutil.Process]: list of child processes, including
                grandchildren.
        """
        return psutil.Process(self.parent_id).children(recursive=True)

    def relentlessly_kill(self, process, timeout=15):
        """
        Tries to kill process, keeps trying until it really is dead.

        Args:
            process (psutil.Process): Process to kill.
            timeout (float): timeout limit in seconds.
        """
        alive, start, pid = True, time.time(), process.pid

        #  log message if terminating parent
        if pid == self.parent_id:
            logging.critical(
                f"It's all over! Terminating parent process with pid {pid}")

        #  keep track of number of kill attempts
        attempts = 0

        #  as long as process remains alive or timeout limit hasn't been reached
        #  attempt to kill it
        while alive:
            try:
                logging_statement = f"Attempting to kill process with pid {pid}"
                if attempts > 0:
                    logging_statement += f". This is attempt number {attempts}."
                logging.critical(logging_statement)

                os.kill(pid, 0)
                process.send_signal(signal.SIGKILL)
            except OSError:
                logging.critical(f"process with pid {pid} killed succesfully")
                alive = False
            finally:
                if timeout and time.time() - start > timeout:
                    alive = False
                    self.timeout = True
                    logging.critical(
                        f"timeout before process with pid {pid} could be"
                        "killed"
                    )
                    if pid == self.parent_id:
                        logging.critical("failed to kill parent process")
                        raise SystemExit(0)
                attempts += 1

    def kill_all(self):
        """
        Kills all processes, including parent and children / granchildren.
        """
        logging.critical("attempting to terminate child processes now")
        for child_process in self.children:
            self.relentlessly_kill(process=child_process)
        logging.critical(
            f"child processes for parent {self.parent_id} terminated")
        parent = psutil.Process(pid=self.parent_id)
        self.relentlessly_kill(process=parent)

def arg_parser():
    """
    Sets arguments for Oligomersoup execute.py CLI.

    Returns:
        object: Namespace object with args as properties.
    """

    #  create ArgumentParser object
    my_parser = argparse.ArgumentParser(
        description="executes a Oligomersoup sequencing workflow",
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

        ripper_file = "unknown"

        #  keep track of which ripper file is being screened
        if kwargs:
            if "ripper_name" in kwargs:
                ripper_file = kwargs["ripper_name"]
            elif "ripper" in kwargs:
                ripper_file = kwargs["ripper"]
        else:
            if args:
                ripper_file = args[0]

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
                f"screening failed after {end} minutes "
                f"with exception: {e}")
        finally:
            if run_info:
                logging.info(
                    f"screen run on {run_info.max_cores} threads with "
                    f"{run_info.init_method} child processes")
        return result

    return log_screen

def exception_handler(timed=False, fatal=False, verbose=False):
    """
    Basic exception handler for simple functions.

    Args:
        func (function): function
        timed (bool): specifies whether to time func execution. If True and
            func is executed successfully, logs time in seconds.
        fatal (bool): specifies whether failure to exit func without Exception
            should be considered effectively as runtime error. If True and
            func call fails, SystemExit is raised
    Returns:
        Any: whatever func returns
    """

    def inner_handler(func):

        @wraps(func)
        def attempt(*args, **kwargs):
            """
            Attempts to execute func.

            Raises:
                SystemExit: raised if fatal == True and func call raises
                    Exception.

            Returns:
                Optional[Any]: Whatever func returns when called. NoneType if
                    func call fails and fatal == False.
            """
            name = func.__name__
            result = None
            if timed:
                start = time.time()
            if verbose:
                logging.info(f"attempting to execute {name}")
            try:
                result = func(*args, **kwargs)
                if timed:
                    logging.info(
                        f"execution of {name} successful after"
                        f"{time.time() - start} seconds.")
                elif verbose:
                    logging.info(f"execution of {name} exited without error")
            except Exception as e:
                logging.exception(f"{name} failed with {e}")
                if fatal:
                    logging.critical(f"critical error upon execution of {name}")
                    raise SystemExit(0)
            return result
        return attempt

    return inner_handler

def check_child_func(size_limit=2E9, pickle_check=False):
    """
    Checks output of functions that are used to pass data from child processes
    to parent. Results of function calls must be less than size_limit and
    pickleable.

    Args:
        size_limit (float, optional): size limit of child function output in
            bytes. Defaults to 2E9.
        pickle_check (bool, optional): specifies whether to check if child
            function output is pickleable.
    Raises:
        SystemExit: raised if result of child function call exceeds size_limit
        SystemExit: raised if result of child function call is not pickleable

    Returns:
        Any or "kill" message: If child process output is less than size_limit
            and pickleable, return whatever function returns. If either of these
            conditions are not met, return string containing "kill" substring.
    """

    def child_checker(child_func):

        @wraps(child_func)
        def execute_child(*args, **kwargs):

            oversized, unpickleable = False, False

            child_func_name = child_func.__name__
            try:
                result = child_func(*args, **kwargs)
            except Exception as e:
                result = -1
                logging.exception(e)
            size = sys.getsizeof(result)

            if size >= size_limit:
                logging.critical(
                    f"Size of output from {child_func_name} in child process is"
                    f"{size}. Size limit is {size_limit} b."
                )
                oversized = True

            #   if pickle_check, make sure output of function is pickleable.
            if pickle_check:
                try:
                    pickle.dumps(result)
                except pickle.PicklingError:
                    logging.critical(
                        f"Output from {child_func_name} is unpickleable")
                    unpickleable = True
            if oversized or unpickleable:
                result = -1
            elif not result:
                return 0
            else:
                print(f"result = {result}")
            return result

        return execute_child
    return child_checker
