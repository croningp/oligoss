"""
This file should be run directly to execute full de novo sequencing
experiments.

"""
import sys
from .extractors.sequence_screening import *
from .insilico.SilicoGenerator import *
from .filehandler import *
from .postprocessing.postprocess import *

def launch_screen(input_parameters_file):
    """
    This function reads an input parameters .json file and decides what
    kind of screening method to use, then performs sequencing screen

    Args:
        input_parameters_file (str): full file path to input parameters .json
            file
    """

    if sys.version_info < (3, 6):
        print("you are running Python version {0}".format(sys.version))
        raise Exception("this package required Python version 3.6 or later")

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

    # load parameters for in silico operations
    silico_params = parameters_dict['silico_parameters']

    # load parameters for data extraction operations
    extractor_params = parameters_dict['extractor_parameters']

    # load important file paths from parameters_dict
    directories = parameters_dict['directories']

    # load location of mass spec data in mzml_ripper .json file format
    ripper_folder = directories['ripper_folder']
    # load output folder for saving data, create if it does not exist
    output_folder = directories['output_folder']
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    write_to_json(
        write_dict=parameters_dict,
        output_json=os.path.join(output_folder, 'run_parameters.json')
    )
    # load parameters for postprocessing operations
    postprocess_parameters = parameters_dict["postprocessing_parameters"]

    # load output folder for saving data, create if it does not exist
    output_folder = directories['output_folder']
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    write_to_json(
        write_dict=parameters_dict,
        output_json=os.path.join(output_folder, 'run_parameters.json')
    )
    # load parameters for postprocessing operations
    postprocess_parameters = parameters_dict["postprocessing_parameters"]

    # load parameters for postprocessing operations
    postprocess_parameters = parameters_dict["postprocessing_parameters"]

    # generate compositional silico dict for screening MS1 EICs
    compositional_silico_dict = generate_MS1_compositional_dict(silico_params)

    if extractor_params["filter"]: 
        
        # filter mass spectra to remove useless data before screening
        pre_filter_rippers(
            ripper_folder=ripper_folder,
            extractor_parameters=extractor_params
        )

        # retrieve file paths to filtered mass spectra
        filtered_ripper_folder = os.path.join(ripper_folder, 'filtered_rippers')

    else:
        filtered_ripper_folder = ripper_folder
        
    filtered_rippers = [
        os.path.join(filtered_ripper_folder, file)
        for file in os.listdir(filtered_ripper_folder)
    ]

    # iterate through filtered ripper .json files
    for filtered_ripper in filtered_rippers:

        # create subfolder for this data set 
        write_folder = os.path.join(
            output_folder, 
            os.path.basename(filtered_ripper)
        )

        if not os.path.exists(write_folder):
            os.mkdir(write_folder)

        # open filtered ripper .json file as dict
        ripper_dict = open_json(filtered_ripper)

        print(f'getting MS1 EICs for {len(compositional_silico_dict)} compositions')
       
        # get MS1 EICs for all compositions
        MS1_EICs = extract_MS1(
            silico_dict=compositional_silico_dict,
            extractor_parameters=extractor_params,
            ripper_dict=ripper_dict
        )

        # save MS1 EICs
        write_MS1_EIC_file(
            input_data_file=filtered_ripper,
            output_folder=write_folder,
            MS1_EICs=MS1_EICs
        )

        # remove compositions that are not present at sufficient abundance
        # from in silico compositionl dict
        compositional_silico_dict = {
            composition: compositional_silico_dict[composition]
            for composition in MS1_EICs
        }
        print(f'compositions = {compositional_silico_dict}')
        
        print(f'generating full MSMS dict for {len(MS1_EICs)} compositions')

        # generate full insilico fragmentation MSMS data for all possible
        # sequences that match compositions found in MS1 EICs
        full_MSMS_silico_dict = generate_MSMS_insilico_from_compositions(
            composition_dict=compositional_silico_dict,
            silico_parameters=silico_params,
            uniques=False
        )
       
        # write full in silico dict to .json file
        write_pre_fragment_screen_sequence_JSON(
            input_data_file=filtered_ripper,
            output_folder=write_folder,
            MSMS_silico_dict=full_MSMS_silico_dict
        )

        print(f'confirming fragments for {len(full_MSMS_silico_dict)} sequences')

        # confirm fragments from silico dict and ripper_dict
        confirmed_fragment_dict = confirm_fragments_sequence_dict(
            silico_dict=full_MSMS_silico_dict,
            ripper_dict=ripper_dict,
            extractor_parameters=extractor_params
        )

        print(f'confirmed fragment dict = {confirmed_fragment_dict}')

        # write confirmed fragment silico dict to .json file
        write_confirmed_fragment_dict(
            input_data_file=filtered_ripper,
            output_folder=write_folder,
            confirmed_fragment_dict=confirmed_fragment_dict
        )

        # assign final confidence scores to all sequences with confirmed
        # fragments
        confidence_scores = assign_confidence_sequences(
            silico_dict=full_MSMS_silico_dict,
            confirmed_dict=confirmed_fragment_dict,
            postprocess_params=postprocess_parameters
        )

        print(f"confidence_scores = {confidence_scores}")
        
        # write confidence scores to .json file
        write_confidence_assignments(
            input_data_file=filtered_ripper,
            output_folder=write_folder,
            confidence_assignments=confidence_scores
        )

        # get minimum confidence for identifying unique fragments from
        # confirmed sequences - SEQUENCES BELOW THIS CONFIDENCE THRESHOLD
        # NOT BE ASSIGNED RETENTION TIMES OR INTENSITIES
        min_confidence_threshold = postprocess_parameters[
            "min_viable_confidence"]

        # generate dict of sequence with high enough confidence scores to
        # assign a retention time and intensity from EICs
        confident_assignments = {
            sequence: {
                "MS1": full_MSMS_silico_dict[sequence]["MS1"],
                "MS2": {
                    frag: masses
                    for frag, masses in full_MSMS_silico_dict[sequence]["MS2"].items()
                    if frag in confirmed_fragment_dict[sequence]
                },
                "peak_list": full_MSMS_silico_dict[sequence]["peak_list"]
            }
            for sequence in confirmed_fragment_dict
            if confidence_scores[sequence] >= min_confidence_threshold
        }

        print(f'{len(confident_assignments)} sequences confirmed with confidence')
        print(f'of {min_confidence_threshold}% or above.')
        print(f'These sequences are now being assigned retention times')

        # find any unique fragments from confirmed fragments of confidently
        # assigned sequences, add these to silico dict of high confidence
        # sequences
        unique_fragment_dict = find_unique_fragments_sequence_dict(
            sequence_dict=confident_assignments
        )

        # write unique fragment dict to output .json file
        write_unique_fragment_dict(
            input_data_file=filtered_ripper,
            output_folder=write_folder,
            unique_fragment_dict=unique_fragment_dict
        )

        # get MS2 EICs for sequences with confirmed unique fragments
        MS2_EICs = extract_MS2(
            silico_dict=unique_fragment_dict,
            extractor_parameters=extractor_params,
            ripper_dict=ripper_dict,
            filter_most_intense=True
        )

        # write MS2 EICs to .json file
        write_MS2_EIC_file(
            input_data_file=filtered_ripper,
            output_folder=write_folder,
            MS2_EICs=MS2_EICs
        )

        # return final retention time and intensity assignments for high
        # confidence sequences
        final_Rt_Is = find_Rt_I_sequences(
            MS1_EICs=MS1_EICs,
            MS2_EICs=MS2_EICs,
            confident_assignments=confident_assignments,
            postprocess_parameters=postprocess_parameters,
            confidence_scores=confidence_scores
        )
        
        # write retention time and intensity assignments to output .csv file
        write_final_retention_time_assignments(
            input_data_file=filtered_ripper,
            output_folder=write_folder,
            final_Rt_I_dict=final_Rt_Is
        )

launch_screen(sys.argv[1])
