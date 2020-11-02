"""
This file contains functions for generating MS1 silico libraries
"""

from .helpers.helpers import (
    ionize_sequence_precursors,
    generate_all_sequences,
    return_monomer_ids_sequence
)


def generate_ms1_ions(
    params,
    polymer,
    sequencing=False
):
    """
    Takes Parameters, Polymer object and yields sequence / composition and MS1
    precursor m/z values in format: (sequence [str], m/z [List[float]]).

    Args:
        params (Parameters): Parameters object.
        polymer (Polymer): Polymer object.
        sequencing (bool, optional): specifies whether to calculare precursors
            for sequences or composition strings; if True, sequence strings are
            yielded (compositions if False). Defaults to False.

    Yields:
        Tuple(str, List[float]): tuple of sequence string and precursor m/z
            values.
    """
    #  generate list of sequences or compositions
    sequences = generate_all_sequences(
        polymer=polymer,
        params=params,
        sequencing=sequencing)

    #  iterate through sequences / compositions and work out MS1 precursor
    #  m/z values, add these to precursor dict
    for seq in sequences:
        if sequencing:
            #  convert sequence to composition for cache
            sorted_seq = "".join(sorted(return_monomer_ids_sequence(
                sequence=seq, return_modified=True, return_set=False)))
        else:
            sorted_seq = seq
        precursors = ionize_sequence_precursors(
            sequence=sorted_seq,
            params=params,
            polymer=polymer)

        yield (seq, precursors)
