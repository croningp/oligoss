"""
OLIGOSS Version 0.0.1

This file is used to run all OLIGOSS experiments.
Artificial Life Team, Cronin Group, 2020
"""
import os
import time
import logging
import faulthandler
import multiprocessing

import concurrent.futures

from logging.config import dictConfig

from .utils.logging_utils.logger_utils import (
    generate_logger_config_dict,
    get_git_info_for_log
)
from .utils.parameter_handlers.parameter_handlers import (
    generate_parameters
)
from .utils.run_utils import (
    arg_parser,
    exception_handler,
    RunInfo,
    log_screening_run,
    check_child_func
)
from .utils.file_io import return_jsons, open_json, FileQueueManager

from .silico.polymer_info.polymer import Polymer
from .silico.ms1_silico import generate_ms1_ions

from .workflows.workflow_helpers.linux_workflows import (
    generate_screen_ms2
)

from .extractors.extractor_classes import RipperDict
from .extractors.filters import prefilter_all, mzml_to_json
from .extractors.data_extraction import extract_MS1_EICs

from .postprocessing.postprocess import postprocess_sequence

#  get Namespace object from CLI
args = arg_parser()

#  enable faulthandler to debug segmentation faults
faulthandler.enable(all_threads=True)

manager = multiprocessing.Manager()

def run_oligoss(input_file, data_folder=None, out_folder=None):
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

@exception_handler(fatal=True, verbose=True)
def exhaustive_screen_multi(params, data_folder, out_folder):
    """
    Exhaustive screen with multiprocessing

    Args:
        params (Parameters): parameters object.
        data_folder (str): ripper data dir.
        out_folder (str): output data dir
    """

    #  retrieve info for OS and number of cores
    run_info = RunInfo(params=params)
    run_info.log_record()

    #  generate and log Polymer object
    polymer = Polymer(params)
    polymer.log_record()

    #  make sure all out files get dumped in "extracted" subfolder
    out_folder = os.path.join(out_folder, "extracted")
    if not os.path.exists(out_folder):
        os.mkdir(out_folder)

    #  convert mzml files to ripper JSONs
    mzml_to_json(
        input_folder=data_folder,
        extractor_parameters=params.extractors
    )

    #  get list of filepaths to ripper JSON files
    rippers = return_jsons(data_folder)

    #  pre - filter ripper JSONs and save output
    filtered_rippers = prefilter_all(
        rippers=rippers,
        extractor_parameters=params.extractors
    )

    if run_info.init_method == "forked":
        screen_rippers_linux(
            filtered_rippers=filtered_rippers,
            params=params,
            polymer=polymer,
            out_folder=out_folder,
            run_info=run_info)

def screen_rippers_linux(
    filtered_rippers,
    params,
    polymer,
    out_folder,
    run_info
):
    """
    Takes a list of filtered rippers, compositional MS1 hits, and performs
    exhaustive screen using Linux-specific functions. NOTE: this version of
    the exhaustive screen assumes fork-safe multiprocessing for reading (not
    writing) ripper objects and does not use MongoDB.

    Args:
        filtered_rippers (List[str]): list of full filepaths to filtered
            ripper files.
        params (Parameters): Parameters object.
        polymer (Polymer): Polymer object.
        out_folder (str): output directory.
        run_info (RunRinfo): RunInfo object.
    """

    #  keep track of how many rippers have been screened
    ripper_count = 0

    #  iterate through the ripper directories
    for ripper in filtered_rippers:

        logging.info("generating MS1 precursors")

        #  generate compositions and their corresponding MS1 m/z values
        compositional_ms1 = generate_ms1_ions(
            params=params,
            polymer=polymer,
            sequencing=False
        )

        ripper_name = os.path.basename(ripper).rstrip(".json")

        #  create ripper output dir
        ripper_output = os.path.join(out_folder, ripper_name)
        if not os.path.exists(ripper_output):
            os.mkdir(ripper_output)

        filequeuer = FileQueueManager(
            output_dir=ripper_output,
            sequence_space=len(polymer)
        )

        print(f"filequeuer ouput = {filequeuer.output_dir}")

        logging.info(
            f"screening ripper file ({ripper}). This is ripper file\n"
            f"{ripper_count + 1} of {len(filtered_rippers)}. Output data for\n"
            f"for this ripper will be saved to {ripper_output}")

        exhaustive_screen_forked(
            ripper=ripper,
            ms1_silico=compositional_ms1,
            params=params,
            polymer=polymer,
            ripper_output=ripper_output,
            run_info=run_info,
            filequeuer=filequeuer
        )
        ripper_count += 1

@log_screening_run
def exhaustive_screen_forked(
    ripper,
    ms1_silico,
    polymer,
    params,
    ripper_output,
    run_info,
    filequeuer
):
    """
    Perform multiprocessed exhaustive screen on single ripper assuming child
    processed are forked (i.e. OS = Linux).

    Args:
        ripper (str): full filepath to ripper.
        polymer (Polymer): Polymer object.
        params (Parameters): Parameters object.
        ripper_output (str): output directory.
        run_info (RunInfo): RunInfo object.
        locks (Dict[str, lock]): string ids and corresponding Lock objects for
            controlling access to output files.
    """
    logging.info(f"opening ripper file: {ripper}")

    #  get Ripper object from ripper dir
    ripper_data = RipperDict(open_json(ripper))

    #  set ripper MS1, MS2 objects for current ripper
    #  NOTE: global keyword is used here to avoid issues with the GLI.
    #  Assigning ripper objects at top-level is fine for single rippers,
    #  but their reassigned values are not picked up by child processes unless
    #  the global keyword is used. I know global is frowned upon but I think
    #  it's justified in this case. Any suggestions for a better solution can
    #  be directed to David Doran or Emma Clarke
    global RIPPER_MS1
    global RIPPER_MS2
    RIPPER_MS1, RIPPER_MS2 = ripper_data.ms1, ripper_data.ms2

    logging.info(
        "launching child processes and screening for presence of precursors"
        "and their isomeric sequences at MS1, then MS2")

    #   multiprocess screening for compositions => screen MS1, generate MS2
    #   silico and then screen MS2 in parallel
    composition_count, problem_count = 0, 0
    with concurrent.futures.ProcessPoolExecutor(
            max_workers=run_info.max_cores) as executor:
        results = [executor.submit(
            full_screen_composition_forked,
            comp, masses, polymer, params, ripper_output, filequeuer)
            for comp, masses in ms1_silico]
        for x in concurrent.futures.as_completed(results):
            composition_count += 1
            if x.result() == 0:
                logging.info(
                    f"composition number {composition_count} screened without "
                    "issue.")
                filequeuer.update()
            else:
                problem_count += 1
                logging.warning(
                    "possible problem with composition number "
                    f"{composition_count}.")
    logging.info(f"Total of {composition_count} compositions screened.")
    if problem_count > 0:
        logging.warning(
            f"Of {composition_count} compositions screened, {problem_count}"
            "potentially have problems.")

    logging.info(filequeuer)

    filequeuer.update()

    logging.info(filequeuer)

@check_child_func()
def full_screen_composition_forked(
    composition,
    precursor_masses,
    polymer,
    params,
    out_folder,
    filequeue
):
    """
    Screens a single ripper for a single composition. Called on as part of
    multiprocessed when screen when child processes have been forked (i.e.
    the default when using Linux).
    """
    #  set out folder for individual composition
    f_basename = f"{os.path.basename(out_folder)}_{composition}"

    logging.info(f"extracting ms1 hit for {composition}")

    #  get MS1 EIC for target composition
    ms1_hit = extract_MS1_EICs(
        MS1_silico={composition: precursor_masses},
        extractor_parameters=params.extractors,
        filename=f_basename,
        MS1_dict=RIPPER_MS1,
        db=False
    )

    #  precursors not present so no point in doing MS2
    if not ms1_hit:
        return None

    #  get maximum intensity for MS1 EIC then remove to free up memory
    max_intensity = ms1_hit[composition][0][1]
    if params.extractors.min_ms1_max_intensity:
        if max_intensity < params.extractors.min_ms1_max_intensity:
            logging.warning(f"no MS1 hit for {composition}")
            return None

    del ms1_hit

    ms2_hits = generate_screen_ms2(
        composition=composition,
        precursors=precursor_masses,
        params=params,
        polymer=polymer,
        out_folder=out_folder,
        ripper_ms2=RIPPER_MS2)

    hits, misses = 0, 0
    for (seq, confirmed, _) in ms2_hits:
        confirmation_data = {
            "sequence": seq,
            "confirmed_fragments": {
                k: v for k, v in confirmed.items() if k != "unconfirmed"},
            "max_intensity": max_intensity,
            "composition": composition,
            "unconfirmed": confirmed["unconfirmed"]
        }
        result = postprocess_sequence(
            sequence_info=confirmation_data, params=params)
        if result:
            hits += 1
            filequeue.queue_data(result)

    if hits == 0:
        logging.info(f"no MS2 hits for composition {composition}")
    else:
        logging.info(
            f"{hits} sequence hits and {misses} misses for composition"
            f"{composition}")


if __name__ == '__main__':
    run_oligoss(
        input_file=args.input,
        data_folder=args.ripper,
        out_folder=args.output
    )
