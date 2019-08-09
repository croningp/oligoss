"""
This file should be run directly to execute full sequence screening
experiments.

"""
from .extractors.sequence_screening import *
from .insilico.SilicoGenerator import *
from .filehandler import *
from .postprocessing.postprocess import *

def launch_screen(input_parameters_file='C:/Users/group/PolymerMassSpec/Examples/InputParams_Test.json'):
    """
    This function reads an input parameters .json file and decides what
    kind of screening method to use, then performs sequencing screen

    Args:
        input_parameters_file (str): full file path to input parameters .json
            file
    """
    # load parameters
    parameters_dict = generate_parameters_dict(input_parameters_file)

    # perform exhaustive screen if that is specified in input parameters
    if parameters_dict['screening_method'] == 'exhaustive':
        exhaustive_screen(parameters_dict)


def exhaustive_screen(
    parameters_dict
):
    """
    This function performs an exhaustive or 'brute force' screen for all
    possible sequences arising from input monomers and constraints set, read
    directly from input parameters .json file

    Args:
        parameters_dict (dict): input parameters dictionary in standard
            input parameters format
    """
    start = time.time()
    # load polymer-specific info - THIS IS CURRENTLY NOT IN USE
    polymer_data = parameters_dict['directories']['polymer_config']

    # load parameters for in silico operations
    silico_params = parameters_dict['silico_parameters']

    # load parameters for data extraction operations
    extractor_params = parameters_dict['extractor_parameters']

    # load important file paths from parameters_dict
    directories = parameters_dict['directories']

    # load location of mass spec data in mzml_ripper .json file format
    ripper_folder = directories['ripper_folder']

    # load parameters for postprocessing operations
    postprocess_parameters = parameters_dict["postprocessing_parameters"]

    # generate compositional silico dict for screening MS1 EICs
    compositional_silico_dict = generate_MS1_compositional_dict(silico_params)

    # filter mass spectra to remove useless data before screening
    pre_filter_rippers(
        ripper_folder=ripper_folder,
        extractor_parameters=extractor_params
    )

    # retrieve file paths to filtered mass spectra
    filtered_ripper_folder = os.path.join(ripper_folder, 'filtered_rippers')
    filtered_rippers = [
        os.path.join(filtered_ripper_folder, file)
        for file in os.listdir(filtered_ripper_folder)
    ]

    # iterate through filtered ripper .json files
    for filtered_ripper in filtered_rippers:

        # open filtered ripper .json file as dict
        ripper_dict = open_json(filtered_ripper)

        print(f'getting MS1 EICs for {len(compositional_silico_dict)} compositions')
        EIC_start = time.time()
        # get MS1 EICs for all compositions
        MS1_EICs = extract_MS1(
            silico_dict=compositional_silico_dict,
            extractor_parameters=extractor_params,
            ripper_dict=ripper_dict
        )
        EIC_duration = time.time()-EIC_start
        print(f'{len(MS1_EICs)} compositions present at sufficient abundance')
        print(f'EIC generator took {EIC_duration/len(compositional_silico_dict)} per composition')

        # remove compositions that are not present at sufficient abundance
        # from in silico compositionl dict
        compositional_silico_dict = {
            composition: compositional_silico_dict[composition]
            for composition in MS1_EICs
        }

        print(f'generating full MSMS dict for {len(MS1_EICs)} compositions')

        # generate full insilico fragmentation MSMS data for all possible
        # sequences that match compositions found in MS1 EICs
        full_MSMS_silico_dict = generate_MSMS_insilico_from_compositions(
            composition_dict=compositional_silico_dict,
            silico_parameters=silico_params,
            uniques=False
        )

        print(f'confirming fragments for {len(full_MSMS_silico_dict)} sequences')

        # confirm fragments from silico dict and ripper_dict
        confirmed_fragment_dict = confirm_fragments_sequence_dict(
            silico_dict=full_MSMS_silico_dict,
            ripper_dict=ripper_dict,
            extractor_parameters=extractor_params
        )
        confirmed_write = f'confirmed_{silico_params["MS1"]["monomers"]}.json'
        write_to_json(confirmed_fragment_dict, os.path.join(ripper_folder, confirmed_write))

        duration = time.time() - start

        print(f'full run for {len(confirmed_fragment_dict)} con seqs = {duration} s')
        # assign final confidence scores to all sequences with confirmed
        # fragments
        confidence_scores = assign_confidence_sequences(
            silico_dict=full_MSMS_silico_dict,
            confirmed_dict=confirmed_fragment_dict,
            postprocess_params=postprocess_parameters
        )

        # get minimum confidence for identifying unique fragments from
        # confirmed sequences - SEQUENCES BELOW THIS CONFIDENCE THRESHOLD
        # NOT BE ASSIGNED RETENTION TIMES OR INTENSITIES
        min_confidence_threshold = postprocess_parameters[
            "min_viable_confidence"]

        # generate dict of sequence with high enough confidence scores to
        # assign a retention time and intensity from EICs
        confident_assignments = {
            sequence: confirmed_fragments
            for sequence, confirmed_fragments in confirmed_fragment_dict.items()
            if confidence_scores[sequence] >= min_confidence_threshold
        }
        print(f'confident_assignments = {confident_assignments}')

        # build trimmed MSMS dict of sequences, fragments and fragment masses
        # for confirmed fragments, excluding signature ions
        confirmed_MSMS_dict = {}

        # iterate through confidently assigned sequences
        for sequence, confirmed_fragments in confident_assignments.items():

            # list in silico signature ions for sequence
            signatures = full_MSMS_silico_dict[sequence]["MS2"]["signatures"]

            # build a trimmed fragment dict for sequence, including only
            # confirmed fragments that are NOT signatures
            confirmed_MSMS_dict[sequence] = {
                frag : full_MSMS_silico_dict[sequence]["MS2"][frag]
                for frag in confirmed_fragments
                if frag not in signatures
            }

        # find unique fragments from those that have been confirmed for each
        # sequence
        unique_confirmed_fragments = find_unique_fragments_sequences_fragment_subset(
            confirmed_MSMS_dict
        )

        # generate MSMS sequence dict for sequences and confirmed UNIQUE
        # fragments for screening and generating MS2 unique EICs for each
        # sequence
        screening_dict = {
            seq: {
                "MS1": full_MSMS_silico_dict[seq],
                "MS2": {
                    frag : masses
                    for frag, masses in confirmed_MSMS_dict[seq]
                    if frag in unique_confirmed_fragments[seq]
                },
                "peak_list": full_MSMS_silico_dict[seq]["peak_list"]
            }
            for seq in unique_confirmed_fragments
        }

        # get MS2 EICs for the most abundant / intense unique fragment
        # for each sequence
        MS2_EICs = extract_MS2(
            silico_dict=screening_dict,
            extractor_parameters=extractor_params,
            ripper_dict=ripper_dict
        )

launch_screen()
