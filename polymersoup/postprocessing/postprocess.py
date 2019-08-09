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
    sequence_EIC_dict,
    postprocess_parameters,
    isobaric_assignment=False
):


    pass
