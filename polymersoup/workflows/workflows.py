"""
This file contains functions for the main Polymersoup experimental workflows
"""
import os
import time
import logging
import concurrent.futures
import multiprocessing

from ..silico.polymer_info.polymer import Polymer
from ..silico.ms1_silico import generate_ms1_ions
from ..silico.silico_handler import (
    get_ms2_silico_dict_from_compositions,
    combine_ms1_ms2_silico_dicts
)
from ..extractors.extractor_classes import RipperDict
from ..extractors.data_extraction import (
    extract_MS1_EICs,
    confirm_all_fragments_concurrent
)
from ..extractors.filters import mzml_to_json, prefilter_all
from ..utils.file_io import (
    write_locked_csv,
    return_jsons,
    open_json,
    write_to_json,
    append_locked_json
)
from ..utils.run_utils import set_max_cores, check_child_process_fork
from ..postprocessing.postprocess_helpers import list_unconfirmed_fragments
from ..postprocessing.postprocess import postprocess_composition

#  get lock object for synchronising access to summary csv, silico and spectral
#   assignment JSON files between processes
CSV_LOCK, JSON_LOCK, SILICO_LOCK = (
    multiprocessing.Lock(),
    multiprocessing.Lock(),
    multiprocessing.Lock()
)

def exhaustive_screen_multi(params, data_folder, out_folder):
    """
    Exhaustive screen with multiprocessing

    Args:
        params (Parameters): parameters object.
        data_folder (str): ripper data dir.
        out_folder (str): output data dir
    """

    #  make sure all out files get dumped in "extracted" subfolder
    out_folder = os.path.join(out_folder, "extracted")
    if not os.path.exists(out_folder):
        print(out_folder)
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

    logging.info("generating Polymer object")

    #  get Polymer object from input parameters
    polymer = Polymer(params_obj=params)

    logging.info("generating MS1 compositional dict")

    #  generate dict of compositions and MS1 precursors for screening
    compositional_ms1_dict = generate_ms1_ions(
        params=params,
        polymer=polymer,
        sequencing=False
    )

    logging.info(
        f"MS1 silico dict generated for {len(compositional_ms1_dict)}\n"
        "compositions".rstrip("\n"))

    #  dump MS1 compositional dict to JSON
    ms1_silico_json = os.path.join(out_folder, "ms1_silico.json")
    write_to_json(
        write_dict=compositional_ms1_dict,
        output_json=ms1_silico_json
    )

    logging.info(f"MS1 silico dict written to {ms1_silico_json}")

    #  perform exhaustive screen on ripper files
    screen_rippers(
        filtered_rippers=filtered_rippers,
        compositional_ms1_dict=compositional_ms1_dict,
        params=params,
        polymer=polymer,
        out_folder=out_folder)


def screen_rippers(
    filtered_rippers: list,
    compositional_ms1_dict: dict,
    params: object,
    polymer: object,
    out_folder: str
):
    """
    Takes a list of filtered rippers, compositional MS1 hits, and performs
    exhaustive sc

    Args:
        filtered_rippers (List[str]): list of full filepaths to filtered
            ripper files.
        compositional_ms1_dict (Dict[str, List[float]]): compositional hits
            and their MS1 precursors.
        params (Parameters): Parameters object.
        polymer (Polymer): Polymer object.
        out_folder (str): output directory.
    """
    #  iterate through the ripper directories
    for ripper in filtered_rippers:

        logging.info(f"screening ripper file ({ripper})")

        #  create ripper output dir
        ripper_output = os.path.join(
            out_folder, os.path.basename(ripper)).replace(".json", "")
        if not os.path.exists(ripper_output):
            os.mkdir(ripper_output)

        logging.info(
            f"output data for screen will be stored in {ripper_output}")

        #  create summary csv for postprocessing results
        ssw = params.postprocess.subsequence_weight
        output_csv = os.path.join(ripper_output, f"{ssw}ssw_summary.csv")
        write_locked_csv(
            fpath=output_csv,
            write_list=[
                "sequence",
                "confidence",
                "confirmed_signatures",
                "max_intensity",
                "composition"],
            lock=CSV_LOCK)
        output_json = os.path.join(ripper_output, "spectral_assignments.json")
        write_to_json(
            write_dict={"source_ripper": ripper},
            output_json=output_json,
            type=str
        )
        write_to_json(
            write_dict={"silico": f"{len(compositional_ms1_dict)}"},
            output_json=os.path.join(ripper_output, "MSMS_silico.json"),
            type=str
        )
        #  get maximum number of cores to use in workflow
        global MAX_CORES
        MAX_CORES = set_max_cores(params)

        #  iterate through rippers, carry out exhaustive screens
        for ripper in filtered_rippers:

            #  check whether OS is Linux or Windows
            if check_child_process_fork():
                exhaustive_screen_forked(
                    ripper=ripper,
                    compositional_ms1_dict=compositional_ms1_dict,
                    polymer=polymer,
                    params=params,
                    ripper_output=ripper_output
                )
            else:
                exhaustive_screen_spawned(
                    ripper=ripper,
                    compositional_ms1_dict=compositional_ms1_dict,
                    params=params,
                    polymer=polymer,
                    ripper_output=ripper_output)

def exhaustive_screen_forked(
    ripper: str,
    compositional_ms1_dict: dict,
    polymer: object,
    params: object,
    ripper_output: str
):
    """
    Perform multiprocessed exhaustive screen on single ripper assuming child
    processed are forked (i.e. Linux).

    Args:
        ripper (str): full filepath to ripper.
        compositional_ms1_dict (Dict[str, List[float]]): compositional hits
            and MS1 precursors.
        polymer (Polymer): Polymer object.
        params (Parameters): Parameters object.
        ripper_output (str): output directory.
    """
    #  get Ripper object from ripper dir
    ripper_data = RipperDict(open_json(ripper))

    #  set ripper MS1, MS2 objects and unique tag for current ripper
    global RIPPER_MS1
    global RIPPER_MS2
    global RIPPER_TAG
    RIPPER_MS1, RIPPER_MS2 = ripper_data.ms1, ripper_data.ms2
    RIPPER_TAG = str(time.time())

    output_csv = os.path.join(
        ripper_output,
        f"{params.postprocess.subsequence_weight}ssw_summary.csv")
    output_json = os.path.join(
        ripper_output,
        "spectral_assignments.json"
    )
    # output_silico_json = os.path.join(
    #     ripper_output,
    #     "MS/MS_silico.json"
    # )

    #   multiprocess screening for compositions => screen MS1, generate MS2
    #   silico and then screen MS2 in parallel
    with concurrent.futures.ProcessPoolExecutor(max_workers=MAX_CORES) \
            as executor:
        results = [executor.submit(
            full_screen_composition_forked,
            comp, masses, polymer, params, ripper_output)
            for comp, masses in compositional_ms1_dict.items()]
        for x in concurrent.futures.as_completed(results):
            if x.result():
                postprocess_composition(
                    hit_info=x.result(),
                    params=params,
                    csv_lock=CSV_LOCK,
                    json_lock=JSON_LOCK,
                    output_csv=output_csv,
                    output_json=output_json)

def exhaustive_screen_spawned(
    ripper: str,
    compositional_ms1_dict: dict,
    params: object,
    polymer: object,
    ripper_output: str
):
    """
    Perform multiprocessed exhaustive screen on single ripper assuming child
    processed are spawned (i.e. Windows).

    Args:
        ripper (str): full filepath to ripper.
        compositional_ms1_dict (Dict[str, List[float]]): compositional hits and
            their associated MS1 precursors.
        params (object): Parameters object.
        polymer (object): Polymer object.
        ripper_output (str): output directory.
    """

    ripper_data = RipperDict(open_json(ripper))

    #  create Manager objects for MS1 and MS2 ripper data for sharing
    #  between processes
    ripper_ms1 = multiprocessing.Manager().dict(ripper_data.ms1)
    ripper_ms2 = multiprocessing.Manager().dict(ripper_data.ms2)
    ripper_tag = os.path.abspath(ripper)

    output_csv = os.path.join(
        ripper_output,
        f"{params.postprocess.subsequence_weight}ssw_summary.csv")

    output_json = os.path.join(
        ripper_output,
        "spectral_assignments.json")

    with concurrent.futures.ProcessPoolExecutor(max_workers=MAX_CORES) \
            as executor:
        results = [executor.submit(
            full_screen_composition_spawned,
            ripper_ms1, ripper_ms2, ripper_tag, comp,
            masses, polymer, params, ripper_output)
            for comp, masses in compositional_ms1_dict.items()]
        for x in concurrent.futures.as_completed(results):
            if x.result():
                postprocess_composition(
                    hit_info=x.result(),
                    params=params,
                    csv_lock=CSV_LOCK,
                    json_lock=JSON_LOCK,
                    output_csv=output_csv,
                    output_json=output_json)


def full_screen_composition_spawned(
    ripper_ms1: object,
    ripper_ms2: object,
    ripper_tag: str,
    composition: str,
    precursor_masses: list,
    polymer: object,
    params: object,
    out_folder: str,
    silico_lock: object
) -> dict:
    """
    Screens a single ripper for a single composition. Called on as part of
    multiprocessed when screen when child processes have been spawned.

    Args:
        ripper_ms1 (object): Manager dict for ripper MS1 data.
        ripper_ms2 (object): Manager dict for ripper MS2 data.
        ripper_tag (str): unique tag for ripper file. Used for caching.
        composition (str): composition string for composition being screened.
        precursor_masses (List[float]): list of m/z values for composition
            precursors.
        polymer (object): Polmyer object.
        params (object): Parameters object.
        out_folder (str): output directory for saving results.
        silico_lock (Lock): Lock object for controlling access to silico JSON.

    Returns:
        (Dict[str, dict], optional): dict of results summarised in following
            format:
            {
                "max_intensity": (float),
                "composition": (str)
                seq: {
                    "confirmed_fragments": (Dict[str, List[float]]),
                    "spectral_assignments": (Dict[str, List[str]])
                }
                for seq in isomeric_sequences
            }
    """
    #  set out folder for individual composition
    f_basename = f"{os.path.basename(out_folder)}_{composition}"

    #  get MS1 EIC for target composition
    ms1_hit = extract_MS1_EICs(
        MS1_silico={composition: precursor_masses},
        extractor_parameters=params.extractors,
        filename=f_basename,
        MS1_dict=ripper_ms1
    )

    #  precursors not present so no point in doing MS2
    if not ms1_hit:
        return None

    ms2_hits = generate_screen_ms2(
        composition=composition,
        precursors=precursor_masses,
        params=params,
        polymer=polymer,
        out_folder=out_folder,
        ripper_ms2=ripper_ms2,
        ripper_tag=ripper_tag,
        silico_lock=SILICO_LOCK)

    if ms2_hits and ms2_hits[0]:
        seq_ms2 = {
            seq: {
                "confirmed_fragments": ms2_hits[0][seq],
                "spectral_assignments": ms2_hits[0][seq]
            }
            for seq in ms2_hits[0]}
        seq_ms2.update({
            "max_intensity": ms1_hit[composition][0][1],
            "composition": composition})

        return seq_ms2

    return None

def full_screen_composition_forked(
    composition: str,
    precursor_masses: list,
    polymer: object,
    params: object,
    out_folder: str
) -> dict:

    """
    Screens a single ripper for a single composition. Called on as part of
    multiprocessed when screen when child processes have been spawned.

    Returns:
        (Dict[str, dict], optional): dict of results summarised in following
            format:
            {
                "max_intensity": (float),
                "composition": (str)
                seq: {
                    "confirmed_fragments": (Dict[str, List[float]]),
                    "spectral_assignments": (Dict[str, List[str]])
                }
                for seq in isomeric_sequences
            }
    """

    #  set out folder for individual composition
    f_basename = f"{os.path.basename(out_folder)}_{composition}"

    #  get MS1 EIC for target composition
    ms1_hit = extract_MS1_EICs(
        MS1_silico={composition: precursor_masses},
        extractor_parameters=params.extractors,
        filename=f_basename,
        MS1_dict=RIPPER_MS1
    )

    #  precursors not present so no point in doing MS2
    if not ms1_hit:
        return None

    ms2_hits = generate_screen_ms2(
        composition=composition,
        precursors=precursor_masses,
        params=params,
        polymer=polymer,
        out_folder=out_folder,
        ripper_ms2=RIPPER_MS2,
        ripper_tag=RIPPER_TAG,
        silico_lock=SILICO_LOCK)

    if ms2_hits and ms2_hits[1]:
        seq_ms2 = {
            seq: {
                "confirmed_fragments": ms2_hits[0][seq],
                "spectral_assignments": ms2_hits[1][seq]
            }
            for seq in ms2_hits[0]}
        seq_ms2.update({
            "max_intensity": ms1_hit[composition][0][1],
            "composition": composition})

        return seq_ms2

    return None

def generate_screen_ms2(
    composition,
    precursors,
    params,
    polymer,
    out_folder,
    ripper_ms2,
    ripper_tag,
    silico_lock
):
    """
    Takes a composition and corresponding MS1 precursors, generates full MS/MS
    silico dict for all sequences isomeric to composition. Screens ripper data
    for all isomeric sequences.

    Args:
        composition (str): composition string.
        precursors (List[float]): precursor m/z values.
        params (Parameters): parameters object.
        polymer (Polymer): polymer object.
        out_folder (str): output data dir.
        ripper_ms2 (dict): ripper MS2 data.
        silico_lock (Lock): Lock object.

    Returns:
        dict, dict: confirmed fragments, spectral matches
    """
    logging.info(
        f"generating isomeric sequences and MS2 fragments for {composition}")
    #  get MS2 silico then full MS/MS silico dict
    silico_dict = get_ms2_silico_dict_from_compositions(
        ms1_hits=[composition],
        params=params,
        polymer=polymer
    )
    silico_dict = combine_ms1_ms2_silico_dicts(
        ms1_silico_dict={composition: precursors},
        ms2_silico_dict=silico_dict
    )

    logging.info(
        f"{len(silico_dict)} sequences to be screened for composition\n"
        f"{composition}".rstrip("\n"))

    #  iterate through isomeric seqs, getting any confirmed fragments and
    #  their spectral matches
    confirmed_fragments, spectral_matches = {}, {}
    for sequence, fragment_dict in silico_dict.items():
        if sequence != "compositions":
            confirmed = confirm_all_fragments_concurrent(
                fragment_dict=fragment_dict,
                ms2_spectra=ripper_ms2,
                params=params,
                ripper_tag=ripper_tag)

            if confirmed and confirmed[0]:
                confirmed_fragments[sequence] = confirmed[0]
                confirmed_fragments[sequence].update(
                    {"unconfirmed": list_unconfirmed_fragments(
                        confirmed_fragments=confirmed[0],
                        silico_dict=silico_dict[sequence]["MS2"]
                    )})

                spectral_matches[sequence] = confirmed[1]
    logging.info(
        f"MS2 fragments screened for sequences isomeric to {composition}")

    #  if sequences have been confimed, save their silico fragments
    if confirmed_fragments:
        append_locked_json(
            fpath=os.path.join(out_folder, "MSMS_silico.json"),
            dump_dict=silico_dict,
            lock=silico_lock)

    return confirmed_fragments, spectral_matches
