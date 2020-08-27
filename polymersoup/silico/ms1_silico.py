"""
This file contains functions for generating MS1 silico libraries
"""

from .helpers.helpers import (
    ionize_sequence_precursors,
    generate_all_sequences,
    return_monomer_ids_sequence
)

def generate_ms1_ions_db(
    params,
    polymer,
    connection,
    sequencing=False,
):
    """
    Generates dict of MS1 precursors from Parameters and Polymers object and
    inserts it into MongoDB database.

    Args:
        params (Parameters): Parameters object
        polymer (Polymer): Polymer object
        sequencing (bool, optional): specifies whether to generate precursor
            dict for monomers or unique compositions. Defaults to False.
        connection (str): mongodb connection string.
    """

    #  generate list of sequences or compositions
    sequences = generate_all_sequences(
        polymer=polymer,
        params=params,
        sequencing=sequencing
    )

    # make sure databse is empty before adding entries (reset from last run)
    ms1db = connection["polymersoup"]["ms1_silico"]

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
        #  get precursor m/z values for compositions and add to precursor dict
        #  for sequence
        ms1db.insert_one({
            "_id": seq,
            "precursors": ionize_sequence_precursors(
                sequence=sorted_seq,
                polymer=polymer,
                params=params)})

    # close mongodb connection
    connection.close()

    #  to save memory, clear cache for memoized function after precursor dict
    #  is complete
    ionize_sequence_precursors.cache_clear()

def generate_ms1_ions(
    params,
    polymer,
    sequencing=False
):

    #  generate list of sequences or compositions
    sequences = generate_all_sequences(
        polymer=polymer,
        params=params,
        sequencing=sequencing)

    #  init dict to store MS1 precursors
    precursors = {}

    #  iterate through sequences / compositions and work out MS1 precursor
    #  m/z values, add these to precursor dict
    for seq in sequences:
        if sequencing:
            #  convert sequence to composition for cache
            sorted_seq = "".join(sorted(return_monomer_ids_sequence(
                sequence=seq, return_modified=True, return_set=False)))
        else:
            sorted_seq = seq
        precursors[seq] = ionize_sequence_precursors(
            sequence=sorted_seq,
            params=params,
            polymer=polymer)

    return precursors
