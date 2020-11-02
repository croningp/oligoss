import logging

from ..utils.run_utils import exception_handler
from .postprocess_helpers import assign_confidence_sequence_concurrent

@exception_handler(fatal=True)
def postprocess_composition(
    hit_info,
    params
):

    #  get composition string
    composition = hit_info["composition"]

    logging.info(f"postprocessing hits for {composition}")

    #  get max intensity @ MS1
    max_intensity = hit_info["max_intensity"]

    #  get subsequence weights for calculating confidence
    ssw = params.postprocess.subsequence_weight

    #  2D list of sequence data for writing to file
    write_array = []

    #  iterate through confirmed sequences, assigning confidence values
    for (seq, info) in hit_info.items():
        if seq not in ["composition", "max_intensity"]:
            fragments = info["confirmed_fragments"]

            confidence = assign_confidence_sequence_concurrent(
                confirmed_fragments=fragments["core"],
                unconfirmed=info["unconfirmed"],
                ssw=ssw,
                params=params
            )

            if confidence > 0:
                core_confirmed = [
                    frag for frag in hit_info[seq][
                        "confirmed_fragments"]["core"]
                    if frag not in [
                        "signatures", "terminal_modifications"]]
                sig_confirmed = [
                    frag for frag in hit_info[seq][
                        "confirmed_fragments"]["signatures"]
                    if frag != "terminal_modifications"]

                write_array.append([
                    seq,
                    confidence,
                    core_confirmed,
                    sig_confirmed,
                    max_intensity,
                    composition])

    logging.info(
        f"about to submit postprocess data for {len(write_array)} sequences"
        "to queue"
    )
    return write_array

def postprocess_sequence(sequence_info, params):

    #  get composition string
    composition = sequence_info["composition"]

    #  get max intensity @ MS1
    max_intensity = sequence_info["max_intensity"]

    #  get subsequence weights for calculating confidence
    ssw = params.postprocess.subsequence_weight

    confidence = assign_confidence_sequence_concurrent(
        confirmed_fragments=sequence_info["confirmed_fragments"]["core"],
        unconfirmed=sequence_info["unconfirmed"],
        ssw=ssw,
        params=params
    )
    if confidence > 0:
        core_confirmed = [
            frag for frag in sequence_info["confirmed_fragments"]
            if frag not in ["signatures", "terminal_modifications"]]
        sig_confirmed = [
            frag for frag in sequence_info["confirmed_fragments"]["signatures"]
            if frag != "terminal_modifications"
        ]
        return [
            sequence_info["sequence"],
            confidence,
            core_confirmed,
            sig_confirmed,
            max_intensity,
            composition
        ]
    return None
