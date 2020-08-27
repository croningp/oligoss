"""
This file contains functions specific to Windows workflows - i.e. workflows
which are only to be run on Windows or Mac. Due to the differences in
multiprocessing between Windows / Mac and Linux it is probably much cleaner to
keep the code separate.

NOTE: at present all Windows / Mac-compatible code can be run on any platform,
but Linux-specific code can only run on Linux.
"""
import os
import sys
import logging
import multiprocessing
import concurrent

from pymongo import MongoClient

#  import utils functions
from ...utils.file_io import (
    write_locked_csv,
    write_to_json,
    open_json,
    mongodb_to_json
)
from ...utils.run_utils import log_screening_run

#  import silico functions
from ...silico.ms1_silico import generate_ms1_ions_db
from ...silico.silico_handler import get_ms2_silico_dict_from_compositions_db

#  import from extractors
from ...extractors.extractor_classes import RipperDict
from ...extractors.data_extraction import (
    extract_MS1_EICs_db,
    confirm_all_fragments_concurrent_db
)

#  import postprocessing functions
from ...postprocessing.postprocess import postprocess_composition

#  retrieve default host for MongoClient - this is platform-dependent
if sys.platform.startswith("linux"):
    MONGO_STR = "mongo"
else:
    MONGO_STR = "localhost"

#  create Lock object for controlling access to output csv file
CSV_LOCK = multiprocessing.Lock()

#  create Lock object for controlling access to output summary JSON file
JSON_LOCK = multiprocessing.Lock()

#  create Lock object for controlling access to output silico JSON file
SILICO_LOCK = multiprocessing.Lock()

def screen_rippers_windows(
    filtered_rippers: list,
    params: object,
    polymer: object,
    out_folder: str,
    run_info: object
):
    """
    Takes a list of filtered rippers, compositional MS1 hits, and performs
    exhaustive screen.
    NOTE: this function is intended for use on Windows and Mac, but it will
    also work on Linux platforms. However, we advise against using this on
    Linux platforms as the Linux-specific screening functions will be much
    faster (typically 2-3x faster).

    Args:
        filtered_rippers (List[str]): list of full filepaths to filtered
            ripper files.
        connection (str): mongodb connection string.
        compositional_ms1_dict (Dict[str, List[float]]): compositional hits
            and their MS1 precursors.
        params (Parameters): Parameters object.
        polymer (Polymer): Polymer object.
        out_folder (str): output directory.
    """

    #  set up MongoDB client
    connection = MongoClient(host=MONGO_STR, port=27017)
    connection.drop_database("polymersoup")

    #  close connection as a precaution - probably unnecessary as we are
    #  unlikely to exceed number of available ports and even if we did, one
    #  extra is unlikely to make the difference. Doesn't hurt though
    connection.close()

    #  generate MS1 precursors and insert into MongoDB
    generate_ms1_ions_db(
        params=params,
        polymer=polymer,
        connection=connection,
        sequencing=False)

    #  dump MS1 compositional collection to JSON
    ms1_silico_json = os.path.join(out_folder, "ms1_silico.json")
    mongodb_to_json(
        collection_name=connection.polymersoup.ms1_silico,
        output_json=ms1_silico_json)

    logging.info(f"MS1 silico dict written to {ms1_silico_json}")

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
                "confirmed_core",
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

        exhaustive_screen_spawned(
            ripper=ripper,
            connection=connection,
            params=params,
            polymer=polymer,
            ripper_output=ripper_output,
            run_info=run_info)

@log_screening_run
def exhaustive_screen_spawned(
    ripper: str,
    connection: str,
    params: object,
    polymer: object,
    ripper_output: str,
    run_info: object
):
    """
    Perform multiprocessed exhaustive screen on single ripper assuming child
    processed are spawned (i.e. Windows).

    Args:
        ripper (str): full filepath to ripper.
        connection (str)): mongodb connection string.
        params (object): Parameters object.
        polymer (object): Polymer object.
        ripper_output (str): output directory.
        run_info (RunInfo): RunInfo object.
    """

    #  open ripper data as RipperDict object
    ripper_data = RipperDict(open_json(ripper))
    ripper_name = os.path.basename(ripper).rstrip(".json")

    #  add ripper ms1 and ripper ms2 data to database
    # has to use for loop so that you can set check_keys to false as decimal
    # place keys, ie m/z values, are invalid
    for k, v in ripper_data.ms1.items():
        connection.polymersoup[f'{ripper_name}_ms1_data'].insert({
            "spectrum_id": k,
            "spectrum": v}, check_keys=False)

    # make sure parent value is a float for future precursor matches
    for k, v in ripper_data.ms2.items():
        connection.polymersoup[f'{ripper_name}_ms2_data'].insert({
            "spectrum_id": k,
            "spectrum": v,
            "parent": float(v["parent"])}, check_keys=False)

    ms1_silico = list(connection["polymersoup"]["ms1_silico"].find())
    connection.close()

    # define filepaths to summary csv and spectral assignments json
    output_csv = os.path.join(
        ripper_output,
        f"{params.postprocess.subsequence_weight}ssw_summary.csv")

    output_json = os.path.join(
        ripper_output,
        "spectral_assignments.json")

    with concurrent.futures.ProcessPoolExecutor(
            max_workers=run_info.max_cores) as executor:
        results = [executor.submit(
            full_screen_composition_spawned, ripper_name, entry["_id"],
            entry["precursors"], polymer, params, ripper_output)
            for entry in ms1_silico]

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
    ripper_name: str,
    composition: str,
    precursor_masses: list,
    polymer: object,
    params: object,
    out_folder: str
) -> dict:
    """
    Screens a single ripper for a single composition. Called on as part of
    multiprocessed when screen when child processes have been spawned.

    Args:
        ripper_name (str): unique tag for ripper file. Used for caching.
        composition (str): composition string for composition being screened.
        precursor_masses (List[float]): list of m/z values for composition
            precursors.
        polymer (object): Polmyer object.
        params (object): Parameters object.
        out_folder (str): output directory for saving results.

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

    connection = MongoClient('localhost', 27017)

    #  set out folder for individual composition
    f_basename = f"{os.path.basename(out_folder)}_{composition}"

    #  get MS1 EIC and maximum intensity for target composition
    ms1_hit, max_int = extract_MS1_EICs_db(
        MS1_silico={composition: precursor_masses},
        extractor_parameters=params.extractors,
        filename=f_basename,
        ripper_name=ripper_name,
        connection=connection)

    #  if precursors not present at MS1, stop the search (no MS2 screening)
    if not ms1_hit:
        return None

    generate_screen_ms2_db(
        composition=composition,
        precursors=precursor_masses,
        params=params,
        polymer=polymer,
        out_folder=out_folder,
        ripper_name=ripper_name,
        connection=connection)

    # get all confirmed fragment entries which match composition
    confirmed_fragments = list(
        connection.polymersoup[
            f"{ripper_name}_ms2_hits"].find({"composition": composition}))

    connection.close()

    if confirmed_fragments:

        # rearrange into postprocessing friendly dict
        seq_ms2 = {
            entry["_id"]: {
                "confirmed_fragments": entry["confirmed_fragments"],
                "unconfirmed": entry["unconfirmed_fragments"],
                "spectral_assignments": entry["spectral_matches"]}
            for entry in confirmed_fragments}
        seq_ms2.update({"max_intensity": max_int, "composition": composition})

        return seq_ms2


def generate_screen_ms2_db(
    composition,
    precursors,
    params,
    polymer,
    out_folder,
    ripper_name,
    connection
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
        ripper_name (str): unique tag for ripper file. Used for caching.
        connection (str): mongodb connection string

    Returns:
        dict, dict: confirmed fragments, spectral matches
    """
    logging.info(
        f"generating isomeric sequences and MS2 fragments for {composition}")
    #  get MS2 silico then full MS/MS silico dict
    get_ms2_silico_dict_from_compositions_db(
        ms1_hits=[composition],
        params=params,
        polymer=polymer,
        connection=connection,
        ripper_name=ripper_name,
        out_folder=out_folder)

    # screen MS2 spectra for fragments corresponding to individual sequences
    # save confirmed fragments and spectral assignments as mongoDB collections
    confirm_all_fragments_concurrent_db(
        composition=composition,
        precursors=precursors,
        params=params,
        connection=connection,
        ripper_name=ripper_name)
