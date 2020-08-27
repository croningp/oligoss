"""
This file contains functions specific to Linux workflows - i.e. workflows
which are only to be run on Linux OS's. Due to the differences in
multiprocessing between Windows / Mac and Linux it is probably much cleaner to
keep the code separate.

NOTE: at present all Windows / Mac-compatible code can be run on any platform,
but Linux-specific code can only run on Linux.
"""
import os
import logging
import multiprocessing
import concurrent.futures

#  import utils functions
from ...utils.file_io import (
    write_locked_csv,
    write_to_json,
    open_json,
    append_locked_json
)
from ...utils.run_utils import log_screening_run

#  import silico functions
from ...silico.ms1_silico import generate_ms1_ions
from ...silico.silico_handler import (
    get_ms2_silico_dict_from_compositions,
    combine_ms1_ms2_silico_dicts
)

#  import from extractors
from ...extractors.extractor_classes import RipperDict
from ...extractors.data_extraction import (
    extract_MS1_EICs,
    confirm_all_fragments_concurrent
)

#  import postprocessing functions
from ...postprocessing.postprocess import postprocess_composition
from ...postprocessing.postprocess_helpers import list_unconfirmed_fragments

#  create Lock object for controlling access to output csv file
CSV_LOCK = multiprocessing.Lock()

#  create Lock object for controlling access to output summary JSON file
JSON_LOCK = multiprocessing.Lock()

#  create Lock object for controlling access to output silico JSON file
SILICO_LOCK = multiprocessing.Lock()

def screen_rippers_linux(
    filtered_rippers: list,
    params: object,
    polymer: object,
    out_folder: str,
    run_info: object
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

    #  generate dict of compositions and corresponding MS1 precursor ions
    compositional_ms1_dict = generate_ms1_ions(
        params=params,
        polymer=polymer,
        sequencing=False
    )

    #  keep track of how many rippers have been screened
    ripper_count = 0

    #  iterate through the ripper directories
    for ripper in filtered_rippers:

        ripper_name = os.path.basename(ripper).rstrip(".json")

        #  create ripper output dir
        ripper_output = os.path.join(out_folder, ripper_name)
        if not os.path.exists(ripper_output):
            os.mkdir(ripper_output)

        logging.info(
            f"screening ripper file ({ripper}). This is ripper file\n"
            f"{ripper_count + 1} of {len(filtered_rippers)}. Output data for\n"
            f"for this ripper will be saved to {ripper_output}")

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

        exhaustive_screen_forked(
            ripper=ripper,
            ms1_silico=compositional_ms1_dict,
            params=params,
            polymer=polymer,
            ripper_output=ripper_output,
            run_info=run_info
        )

@log_screening_run
def exhaustive_screen_forked(
    ripper: str,
    ms1_silico: dict,
    polymer: object,
    params: object,
    ripper_output: str,
    run_info: object
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
    """

    #  get Ripper object from ripper dir
    ripper_data = RipperDict(open_json(ripper))

    #  set ripper MS1, MS2 objects and unique tag for current ripper
    #  NOTE: global keyword is used here to avoid issues with the GLI
    #  assigning ripper objects at top-level is fine for single rippers,
    #  but their reassigned values are not picked up by child processes unless
    #  the global keyword is used.
    global RIPPER_MS1
    global RIPPER_MS2
    RIPPER_MS1, RIPPER_MS2 = ripper_data.ms1, ripper_data.ms2

    #  write summary files for confidence csv and spectral assignment JSON
    output_csv = os.path.join(
        ripper_output,
        f"{params.postprocess.subsequence_weight}ssw_summary.csv")
    output_json = os.path.join(
        ripper_output,
        "spectral_assignments.json")

    #  write silico JSON
    write_to_json(
        ripper,
        os.path.join(ripper_output, "MSMS_silico.json"),
        type=str)

    #   multiprocess screening for compositions => screen MS1, generate MS2
    #   silico and then screen MS2 in parallel
    with concurrent.futures.ProcessPoolExecutor(
            max_workers=run_info.max_cores) as executor:
        results = [executor.submit(
            full_screen_composition_forked,
            comp, masses, polymer, params, ripper_output)
            for comp, masses in ms1_silico.items()]
        for x in concurrent.futures.as_completed(results):
            if x.result():
                postprocess_composition(
                    hit_info=x.result(),
                    params=params,
                    csv_lock=CSV_LOCK,
                    json_lock=JSON_LOCK,
                    output_csv=output_csv,
                    output_json=output_json)


def generate_screen_ms2(
    composition,
    precursors,
    params,
    polymer,
    out_folder,
    ripper_ms2
):
    """
    Takes a composition and corresponding MS1 precursors, generates full MS/MS
    silico dict for all sequences isomeric to composition. Screens ripper data
    for all isomeric sequences.

    Args:
        composition (str): target composition.
        precursors (List[float]): list of precursor m/z values.
        params (Parameters): Parameters object.
        polymer (Polymer): Polymer object.
        out_folder (str): output folder for saving data.
        ripper_ms2 (Dict[str, dict]): ripper MS2 dict.

    Returns:
        dict, dict : confirmed_fragments, spectral_matches
    """

    logging.info(
        f"generating isomeric sequences and MS2 fragments for {composition}")
    #  get MS2 silico then full MS/MS silico dict
    silico_dict = get_ms2_silico_dict_from_compositions(
        ms1_hits=[composition],
        params=params,
        polymer=polymer)
    silico_dict = combine_ms1_ms2_silico_dicts(
        ms1_silico_dict={composition: precursors},
        ms2_silico_dict=silico_dict)

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
                params=params)

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
            lock=SILICO_LOCK)

    return confirmed_fragments, spectral_matches


def full_screen_composition_forked(
    composition: str,
    precursor_masses: list,
    polymer: object,
    params: object,
    out_folder: str
) -> dict:
    """
    Screens a single ripper for a single composition. Called on as part of
    multiprocessed when screen when child processes have been forked (i.e.
    the default when using Linux).

    Returns:
        Optional[Dict[str, dict]]: dict of results summarised in following
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
        NOTE: if no hits are found at MS2, NoneType is returned.
    """

    #  set out folder for individual composition
    f_basename = f"{os.path.basename(out_folder)}_{composition}"

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

    ms2_hits = generate_screen_ms2(
        composition=composition,
        precursors=precursor_masses,
        params=params,
        polymer=polymer,
        out_folder=out_folder,
        ripper_ms2=RIPPER_MS2)

    if ms2_hits and ms2_hits[1]:
        seq_ms2 = {
            seq: {
                "confirmed_fragments": {
                    k: v for k, v in ms2_hits[0][seq].items()
                    if k != "unconfirmed"
                },
                "spectral_assignments": ms2_hits[1][seq],
                "unconfirmed": ms2_hits[0][seq]["unconfirmed"]
            }
            for seq in ms2_hits[0]}
        seq_ms2.update({
            "max_intensity": ms1_hit[composition][0][1],
            "composition": composition})

        return seq_ms2

    return None
