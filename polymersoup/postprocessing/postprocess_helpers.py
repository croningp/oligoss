"""
This file contains core operations required by postprocessing module
"""

def get_subsequence_coverage(
    sequence_length,
    confirmed_fragments,
    core_fragment_series,
    optional_core_fragments,
    excluded_core_fragments
):
    """
    This functions calculates the continuous subsequence coverage of a sequence
    hit for each core_fragment_series and returns value as % continuous
    coverage

    Args:
        sequence_length (int): total sequence length (in monomer units)
        confirmed_fragments (list): list of confirmed fragment ids
        core_fragment_series (list): list of fragment one letter codes for
                        series that are to be used in calculation
        optional_core_fragments (list): list of specific fragment ids to
                        exclude from calculation ONLY IF they are not in the
                        list of confirmed fragments. Include in the calculation
                        if they are in the list
        excluded_core_fragments (list): list of specific fragment ids to
                        exclude from subsequence coverage in ALL cases
    """

    pass
