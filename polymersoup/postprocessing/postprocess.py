"""
This file contains functions to pass on postprocessing parameters and
execute postprocessing by calling on core functions from postprocess_helpers.py

"""
from .postprocess_helpers import *


def assign_confidence_sequences(
    silico_dict,
    confirmed_dict,
    postprocess_params
):
    # remove sequences from in silico dict that do not have any confirmed
    # fragments
    silico_dict = {
        sequence : subdict
        for sequence, subdict in silico_dict.items()
        if sequence in confirmed_dict
    }

    # initiate dict to store confidence scores for sequences
    confidence_dict = {}

    # iterate through sequences that have confirmed fragments, assigning each
    # a confidence score
    for sequence, confirmed_fragments in confirmed_dict.items():

        # assign sequence a confidence score, using constraints specified in
        # input postprocess parameters from input parameters .json file
        confidence_assignment = assign_confidence_sequence(
            insilico_fragments=silico_dict[sequence]["MS2"].keys(),
            confirmed_fragments=confirmed_fragments,
            core_fragment_series=postprocess_params["core_linear_series"],
            optional_fragments=postprocess_params["optional_core_frags"],
            exclude_fragments=postprocess_params["excluded_fragments"],
            essential_fragments=postprocess_params["essential_fragments"],
            essential_fragment_cap=postprocess_params[
                "dominant_signature_cap"],
            sequence_coverage_weight=postprocess_params["subsequence_weight"]
        )

        # add sequence and its final confidence score to confidence_dict
        confidence_dict[sequence] = confidence_assignment

    # return dictionary of sequences and their confidence scores
    return confidence_dict

def find_Rt_I_sequences(
    MS1_EICs,
    MS2_EICs,
    confident_assignments,
    postprocess_parameters,
    confidence_scores,
    return_confidence=True
):
    """
    Takes dicts of MS1 and MS2 EICs, confidently assigned sequences, and fits
    a retention time and intensity to each confidently assigned sequence

    Args:
        MS1_EICs (dict): dictionary of COMPOSITIONS and corresponding MS1 EICs
        MS2_EICs (dict): dictionary of SEQUENCES and corresponding MS2 EICs
        confident_assignments (dict): dictionary of sequences and their
            in silico MS1 and MS2 data for sequences that have been confirmed
            with high confidence
        postprocess_parameters (dict): dictionary of postprocess parameters
            passed on from input parameters .json file
        confidence_scores (dict): dictionary of confirmed sequences and 
            assigned confidence scores 
    Kewyword Args:
        return_confidence (bool): specifies whether to return sequences' 
            confidence score in output. Defaults to True 

    Raises:
        Exception: [description]

    Returns:
        [type]: [description]
    """
    compositions = {}

    for sequence in confident_assignments:
        composition = "".join(sorted(sequence))

        if composition not in compositions:
            compositions[composition] = [sequence]

        else:
            compositions[composition].append(sequence)

        if composition not in MS1_EICs:
            raise Exception(f'something has gone wrong with {sequence}')


    final_Rt_Is = {}

    for seq in MS2_EICs:
        seq_Rt_I = get_Rt_I_from_ms2_EIC(
            MS1_EIC = MS1_EICs["".join(sorted(seq))],
            MS2_EIC=MS2_EICs[seq],
            ms2_Rt_bin=postprocess_parameters["ms2_Rt_bin"],
            flexible_ms2_rt=postprocess_parameters["ms2_Rt_flexible"]
        )
        final_Rt_Is[seq] = seq_Rt_I

    for composition in compositions:

        sequences = [
            seq for seq in compositions[composition]
            if (seq not in MS2_EICs
                and seq in confident_assignments)
        ]


        seq_Rt_Is = get_Rt_I_from_ms1_EIC(
            EIC=MS1_EICs[composition],
            Rt_bin=postprocess_parameters["Rt_bin"],
            backup_Rt_bin=postprocess_parameters["backup_Rt_bin"],
            n_targets=len(sequences),
            min_relative_intensity=postprocess_parameters[
                "min_relative_intensity"]
        )

        print(f'following data points found for {sequences}: {seq_Rt_Is}')
        for i in range(0, len(seq_Rt_Is)):
            final_Rt_Is[sequences[i]] = seq_Rt_Is[i]

        final_Rt_Is.update(
            {
                seq: ['unassigned', 'unassigned']
                for seq in sequences
                if seq not in final_Rt_Is
            }
        )

    # check whether to return sequence confidence scores; if so, add to output
    if return_confidence:
        for seq in final_Rt_Is:
            final_Rt_Is[seq].append(
                confidence_scores[seq]
            )

    return final_Rt_Is
