"""
This file contains functions for passing on silico parameters to insilico
functions for generating theoretical MS and MSMS sequence data

"""

from .MS1_silico import *
from .MS2_silico import *
from .silico_helpers.insilico_helpers import *

def generate_MSMS_sequence_dict(
    parameters_file="C:/Users/group/PolymerMassSpec/Examples/InputParams_Test.json",
    write=True,
    output_folder=None,
    ms1=True,
    ms2=True
):
    """
    This function is used to generate a full MS and MSMS in silico sequence
    dict in format: {
                        seq: {
                            'MS1' : [m/z...],
                            'MS2': {
                                'frag': [m/z...],
                                'signatures': {
                                    'signature': [m/z...]
                                },
                                'unique_fragments' ['unique_frag'...]
                            }
                        }
        where seq = sequence string, m/z = m/z value of a given ion, 'frag' =
        MS2 fragment id string, 'signature' = MS2 fragment id string for a
        monomer-specific signature fragment, 'unique_frag' = MS2 fragment id
        string for a fragment with m/z unique to that sequence (i.e. not found
        in MS2 fragments of other, isobaric sequences)

    Args:
        parameters_file (str, optional): full file path to input parameters
            .json file. Defaults to "C:/Users/group/PolymerMassSpec/Examples/InputParams_Test.json".
        write (bool, optional): specifies whether to write sequence dict to
            file. Defaults to True.
        output_folder (str, optional): full file path to output folder to put
            output file; if None, file is dumped in same directory as
            parameters_file. Defaults to None.
        ms1 (bool, optional): specifies whether to generate theoretical MS1
            data for sequences. Defaults to True.
        ms2 (bool, optional): specifies whether to generate theoretical MS2
            data for sequences. Defaults to True.

    Returns:
        full_MSMS_dict (dict): MSMS sequence dict, the format of which
            is described in detail above
    """
    # reads in silico parameters from input parameters json
    silico_params = return_silico_parameters(parameters_file)

    if ms1:
        start = time.time()
        ms1 = silico_params["MS1"]

        MS1 = generate_ms1_mass_dictionary_adducts_losses(
            ms1["monomers"],
            ms1["max_length"],
            ms1["ms1_adducts"],
            silico_params["mode"],
            ms1["min_z"],
            ms1["max_z"],
            ms1["losses"],
            ms1["max_neutral_losses"],
            ms1["loss_products_adducts"],
            ms1["min_length"],
            ms1["chain_terminators"],
            ms1["universal_rxn"],
            ms1["terminal_tags"]["0"],
            ms1["terminal_tags"]["-1"]
        )
        end = time.time()
        time_elapsed = end-start
        print(f'for {len(MS1)} sequences, time taken for MS1 = {time_elapsed} seconds')
        print(f'time taken per sequence for MS1 = {time_elapsed/len(MS1)} seconds')

    start = time.time()
    ms2 = silico_params["MS2"]
    MS2 = generate_ms2_mass_dictionary(
            MS1.keys(),
            ms2["fragment_series"],
            ms2["ms2_adducts"],
            silico_params["mode"],
            ms2["add_signatures"],
            ms2["signatures"],
            ms2["ms2_losses"],
            ms2["ms2_max_neutral_losses"],
            ms2["ms2_loss_products_adducts"],
            ms2["uniques"]
    )
    end = time.time()
    time_elapsed = end-start
    print(f'for {len(MS2)} sequences, time taken for MS2 = {time_elapsed} seconds')
    print(f'time taken per sequence for MS2 = {time_elapsed/len(MS2)} seconds')
    full_MSMS_dict = {}

    for sequence in MS1:
        full_MSMS_dict[sequence] = {
            'MS1': MS1[sequence],
            'MS2': MS2[sequence]
        }

    full_MSMS_dict = add_peak_lists_massdict(full_MSMS_dict)

    if write:
        if not output_folder:
            output_folder = os.path.dirname(parameters_file)

        write_file = generate_insilico_writefile_string(
            output_folder,
            silico_params
        )

        write_to_json(full_MSMS_dict, write_file)

    return full_MSMS_dict


def generate_MS1_compositional_dict(
    silico_parameters
):
    """
    This function generates a dictionary of compositions (not sequences) and
    associated MS1 ions for all possible compositions that can arise from
    constraints specified in input parameters .json file

    Args:
        silico_parameters (dict): dictionary of in silico parameters passed on
            from input parameters .json file

    Returns:
        MS1_dict (dict): dictionary of compositions and their corresponding
            MS1 m/z values in format : {composition: [m/z...]...}
    """

    # load MS1 in silico parameters
    ms1_params = silico_parameters["MS1"]

    # pass on MS1 silico parameters to MS1_Silico.py function to create the
    # composition MS1 dictionary
    MS1_dict = generate_ms1_mass_dictionary_adducts_losses(
        monomers=ms1_params["monomers"],
        max_length=ms1_params["max_length"],
        adducts=ms1_params["ms1_adducts"],
        mode=silico_parameters["mode"],
        min_z=ms1_params["min_z"],
        max_z=ms1_params["max_z"],
        losses=ms1_params["losses"],
        max_total_losses=ms1_params["max_neutral_losses"],
        loss_product_adducts=ms1_params["loss_products_adducts"],
        min_length=ms1_params["min_length"],
        chain_terminators=ms1_params["chain_terminators"],
        start_tags=ms1_params["terminal_tags"]["0"],
        end_tags=ms1_params["terminal_tags"]["-1"],
        sequencing=False,
        isobaric_targets=ms1_params["isobaric_targets"]
    )

    return MS1_dict

def generate_MSMS_insilico_from_compositions(
    composition_dict,
    silico_parameters,
    uniques=False
):
    """
    This function takes a dictionary of product compositions and associated
    MS1 ions, and returns a full_MSMS_sequence dictionary containing in silico
    data for both MS1 ions and MS2 fragments for every possible sequence that
    matches a composition in the input composition_dict

    Args:
        composition_dict (dict): dictionary of compositions and corresponding
            MS1 ion m/z values
        silico_parameters (dict): dictionary of in silico parameters passed on
            from input parameters .json file
        uniques (bool, optional): specifies whether to find unique MS2
            fragmentsfor each sequence and return these in the final sequence
            dict.
            Defaults to False.

    Returns:
        [type]: [description]
    """
    full_MSMS_sequence_dict = {}

    ms1_params = silico_parameters["MS1"]
    ms2_params = silico_parameters["MS2"]


    for composition, MS1_ions in composition_dict.items():
        monomers = list(set([c for i, c in enumerate(composition)]))
        sequence_length = len(composition)

        start_tags = ms1_params["terminal_tags"]["0"]
        end_tags = ms1_params["terminal_tags"]["-1"]

        if start_tags:
            sequence_length -= 1
        if end_tags:
            sequence_length -= 1

        sequences = generate_all_sequences(
            monomers,
            max_length=sequence_length,
            min_length=sequence_length,
            sequencing=True,
            chain_terminators=ms1_params["chain_terminators"],
            start_tags=start_tags,
            end_tags=end_tags
        )

        sequences = list(filter(
            lambda seq: sorted(seq) == sorted(composition),
            sequences
        ))

        ms2_dict = generate_ms2_mass_dictionary(
            sequences=sequences,
            fragment_series=ms2_params["fragment_series"],
            adducts=ms2_params["ms2_adducts"],
            mode=silico_parameters["mode"],
            add_signatures=ms2_params["add_signatures"],
            signatures=ms2_params["signatures"],
            losses=ms2_params["ms2_losses"],
            max_total_losses=ms2_params["ms2_max_neutral_losses"],
            loss_product_adducts=ms2_params["ms2_loss_products_adducts"],
            uniques=False
        )

        full_MSMS_sequence_dict.update(
            {sequence: {
                "MS1": MS1_ions,
                "MS2": ms2_dict[sequence]
            }
            for sequence in ms2_dict
        })

        full_MSMS_sequence_dict = add_peak_lists_massdict(
            massdict=full_MSMS_sequence_dict
        )

    return full_MSMS_sequence_dict

def find_unique_fragments_sequence_dict(
    sequence_dict
):
    """
    Takes a dictionary of sequences and associated MS2 fragments, and returns
    a dictionary of sequences and lists of fragments that are unique to each
    sequences. This function is to be used to identify unique fragments from
    subsets of fragments, not full theoretical MSMS in silico sequence dicts.
    Intended use is for finding unique fragments from fragments that have been
    confirmed as being associated with sequences by functions in extractors
    and postprocessing

    Args:
        sequence_dict (dict): dictionary of sequences and MS2 fragments in
            format: {
                seq: {
                    'frag': [m/z...]
                }
                ...
            }
        where seq = sequence, 'frag' = fragment id, and m/z = m/z value of
        fragment ion(s)

    Returns:
        unique_fragment_dict (dict): dictionary of sequences and associated
            fragment ids for fragments that are unique to each sequence
    """

    # group sequences by composition for comparison of isobaric sequences
    isobaric_sequence_dict = generate_dict_isobaric_sequences(
        sequences=sequence_dict.keys()
    )

    # initiate dict to store sequences and their unique fragments for output
    unique_fragment_dict = sequence_dict

    # iterate through isobaric sequence dict to access lists of isobaric
    # sequences, which have the same composition
    for isobaric_sequences in isobaric_sequence_dict.values():

        # if there is only one sequence of this composition, there is no
        # need to look for unique fragments
        if len(isobaric_sequences) == 1:
            unique_fragment_dict[isobaric_sequences[0]].update(
                {"unique_fragments": []}
            )

        else:
            # generate mini ms2 fragment dict for current set of isobaric
            # sequences
            ms2_dict = {
                seq: sequence_dict[seq]["MS2"]
                for seq in isobaric_sequences
            }

            # find unique fragments for each sequence in ms2_dict
            unique_fragment_isobaric_dict = find_unique_fragments_isobaric_set(
                ms2_dict
            )

            for seq in ms2_dict:
                unique_fragment_dict[seq].update(
                    {"unique_fragments": unique_fragment_isobaric_dict[seq]}
                )
    return unique_fragment_dict
