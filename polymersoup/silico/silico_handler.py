"""
This file contains silico functions to be called directly in sequencing
workflows.
"""
from .ms2_silico import (
    add_signature_ions_sequence,
    build_fragment_dict_for_sequences,
    build_single_fragment
)
from .helpers.helpers import get_isomeric_seqs

def get_ms2_silico_dict_from_compositions(
    ms1_hits,
    params,
    polymer
):
    """
    Takes a list of compositions that are present at MS1, returns full MS2
    silico dict for all sequences that could match one or more composition.

    Args:
        ms1_hits (List[str]): list of composition strings
        params (Parameters): Parameters object
        polymer (Polymer): Polymer object

    Returns:
        Dict[str, dict]: dict of sequences and corresponding MS2 fragments
    """

    #  init dict to store MS2 fragment silico ions for
    #  all sequences matching compositions in ms1_hits
    ms2_silico_dict = {}

    for composition in list(filter(lambda x: x != "retention_times", ms1_hits)):

        #  add signatures here for compositions
        signature_ions = add_signature_ions_sequence(
            sequence=composition,
            params=params,
            polymer=polymer
        )

        #  get all sequences that match target composition
        sequences = get_isomeric_seqs([composition])

        #  for linear fragment series, generate dict MS/MS dict of all sequences
        for linear_series in polymer.fragment_info:

            #  get fragment dict for all sequences
            series_dict = build_fragment_dict_for_sequences(
                sequences=sequences,
                params=params,
                fragment_series=linear_series,
                polymer=polymer
            )

            #  update ms2_silico_dict for each sequence and its fragments
            for sequence in sequences:
                if sequence not in ms2_silico_dict:

                    ms2_silico_dict[sequence] = series_dict[sequence]

                    ms2_silico_dict[sequence].update({
                        "signatures": signature_ions,
                        "composition": composition
                    })
                else:
                    ms2_silico_dict[sequence].update(series_dict[sequence])

            #  to save memory, clear cache of fragment builder function
            build_single_fragment.cache_clear()

    return ms2_silico_dict

def combine_ms1_ms2_silico_dicts(
    ms1_silico_dict,
    ms2_silico_dict
):
    """
    Takes an MS1 silico dict, an MS2 silico dict and combines them into a full
    MS/MS silico dict.

    Args:
        ms1_silico_dict (Dict[str, List[float]]): dict of MS1 composition
            strings and precursor ion m/z values
        ms2_silico_dict (Dict[str, dict]): dict of sequence strings and
            corresponding MS2 fragment ions

    Returns:
        Dict[Tuple[str], dict]: full_msms_silico dict in format: {
            ("sequence", "composition"): {
                "MS1": [m/z...],
                "MS2": {
                    "frag": [m/z...],
                    ...
                    "signatures": {
                        "signature_frag": [m/z...],
                        ...
                    }
                },
                "peak_list: [m/z...]
            }
        }
    """

    #  get full MS/MS silico dict of precursors and MS2 fragments
    full_msms_silico = {}

    for seq, seq_info in ms2_silico_dict.items():

        full_msms_silico[seq] = {
            "MS1": ms1_silico_dict[seq_info["composition"]],
            "MS2": seq_info
        }

    #  iterate through all sequences in silico dict, adding a peak list for each
    #  This is a list of ALL unique ions (both MS1 and MS2) for a sequence
    for seq in full_msms_silico:
        peak_list = [peak for peak in full_msms_silico[seq]["MS1"]]
        for fragment, fragment_masses in full_msms_silico[seq]["MS2"].items():
            if fragment not in ["signatures", "composition"]:
                peak_list.extend(fragment_masses)
            elif fragment == "signatures":
                if full_msms_silico[seq]["MS2"]["signatures"]:
                    for signature_ions in full_msms_silico[seq][
                            "MS2"]["signatures"].values():
                        peak_list.extend(signature_ions)
        full_msms_silico[seq]["peak_list"] = sorted(list(set(peak_list)))

    #  update msms silico dict for list of composition strings
    full_msms_silico.update({
        "compositions": [key for key in ms1_silico_dict]
    })

    return full_msms_silico
