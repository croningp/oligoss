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
    Generates dict of MS1 precursors from Parameters and Polymers object.

    Args:
        params (Parameters): Parameters object
        polymer (Polymer): Polymer object
        sequencing (bool, optional): specifies whether to generate precursor
            dict for monomers or unique compositions. Defaults to False.

    Returns:
        Dict[str, List[float]]: dict of sequence strings and associated list of
            precursor m/z values.
    """
    #  generate list of sequences or compositions
    sequences = generate_all_sequences(
        polymer=polymer,
        params=params,
        sequencing=sequencing
    )

    #  init dict to store ms1 precursor m/z values
    ms1_precursor_dict = {}

    #  iterate through sequences, calculating precursors or retrieving from
    #  cache via lru_cache decorated function
    for seq in sequences:

        if sequencing:
            #  convert sequence to composition for cache
            sorted_seq = "".join(sorted(return_monomer_ids_sequence(
                sequence=seq, return_modified=True, return_set=False)))
        else:
            sorted_seq = seq

        #  get precursor m/z values for compositions and add to precursor dict
        #  for sequence
        ms1_precursor_dict[seq] = ionize_sequence_precursors(
            sequence=sorted_seq,
            polymer=polymer,
            params=params
        )

    #  to save memory, clear cache for memoized function after precursor dict
    #  is complete
    ionize_sequence_precursors.cache_clear()

    return ms1_precursor_dict
