"""
This file contains functions specific to Linux workflows - i.e. workflows
which are only to be run on Linux OS's. Due to the differences in
multiprocessing between Windows / Mac and Linux it is probably much cleaner to
keep the code separate.

NOTE: at present all Windows / Mac-compatible code can be run on any platform,
but Linux-specific code can only run on Linux.
"""
import logging

#  import utils functions
from ...utils.run_utils import exception_handler

#  import silico functions
from ...silico.silico_handler import generate_msms_dicts_for_composition

#  import from extractors
from ...extractors.data_extraction import confirm_all_fragments_concurrent

#  import postprocessing functions
from ...postprocessing.postprocess_helpers import list_unconfirmed_fragments

@exception_handler(fatal=True, timed=True)
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

    Yields:
        str, dict, dict: sequence (str), confirmed fragments (dict),
            spectral matches (dict)
    """

    logging.info(
        f"generating isomeric sequences and MS2 fragments for {composition}")

    #  get MS2 silico then full MS/MS silico dict
    silico_dict = generate_msms_dicts_for_composition(
        composition=composition,
        precursor_ions=precursors,
        polymer=polymer,
        params=params)

    #  iterate through isomeric seqs, getting any confirmed fragments and
    #  their spectral matches
    seq_count, confirmed_count = 0, 0

    for (sequence, fragment_dict) in silico_dict:
        seq_count += 1
        if sequence != "compositions":
            confirmed = confirm_all_fragments_concurrent(
                fragment_dict=fragment_dict,
                ms2_spectra=ripper_ms2,
                params=params)

            if confirmed[0]:
                if confirmed[0]["core"] or confirmed[0]["signatures"]:
                    confirmed_count += 1
                    confirmed_fragments = confirmed[0]
                    unconfirmed = list_unconfirmed_fragments(
                        confirmed_fragments=confirmed[0],
                        silico_dict=fragment_dict["MS2"])
                    spectral_matches = confirmed[1]
                    confirmed_fragments.update(
                        {"unconfirmed": unconfirmed})

                    yield sequence, confirmed_fragments, spectral_matches
    logging.info(
        f"of {seq_count} sequences for composition {composition}, "
        f"{confirmed_count} sequences were confirmed at MS2"
    )
