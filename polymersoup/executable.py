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
        print(f"you are running Python version {sys.version}")
        raise Exception("this package requires Python version 3.6 or later")

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

    # generate compositional silico dict for screening MS1 EICs
    compositional_silico_dict = generate_MS1_compositional_dict(silico_params)

    # apply any pre-filtering steps to ripper data 
    filtered_rippers = apply_filters(
        extractor_parameters=extractor_params,
        ripper_folder=ripper_folder
    )

    # if data_extraction set to True, extract data 
    if parameters_dict['data_extraction']:

        extracted_data_dirs = standard_extraction(
            rippers=filtered_rippers,
            output_folder=output_folder,
            silico_parameters=silico_params,
            compositional_silico_dict=compositional_silico_dict,
            extractor_parameters=extractor_params
        )

    else:

        extracted_data_dirs = [
            os.path.join(output_folder, subdir) 
            for subdir in os.listdir(output_folder)
            if not (os.path.isfile(subdir) 
            and subdir.find('run_parameters') == -1)
        ]

        extracted_data_dirs = [
            os.path.join(data_dir, 'extracted')
            for data_dir in extracted_data_dirs
        ]

        print(f'extracted data dirs = {extracted_data_dirs}')
    
    if parameters_dict['postprocess']:

        for extracted_data in extracted_data_dirs:

            standard_postprocess(
                extracted_data_folder=extracted_data,
                postprocess_parameters=postprocess_parameters
            )

def apply_filters(
    extractor_parameters,
    ripper_folder
):
    """
    Applies any pre-run filters to ripper data and returns list of full 
    file paths to ripper jsons
    
    Args:
        extractor_parameters (dict): dictionary of extractor parameters, 
            retrieved from input parameters JSON file
        ripper_folder (str): full file path to folder containing mzml ripper
            JSON files
    
    Returns:
        list: list of full file paths to mzml ripper JSON files 
    """
    
    # if filter set to False, return list of full file paths to unfiltered
    # rippers 
    if not extractor_parameters["pre_run_filter"]:

        return [
            os.path.join(ripper_folder, ripper) 
            for ripper in os.listdir(ripper_folder)
            if ripper.endswith('.json')
        ]
    
    # if filter set to True, apply filter and return list of full file paths
    # to newly created filtered JSON files 
    return pre_filter_rippers(
        ripper_folder=ripper_folder,
        extractor_parameters=extractor_parameters
        )


def standard_extraction(
    rippers,
    extractor_parameters,
    silico_parameters,
    compositional_silico_dict,
    output_folder

):

    # if rippers data folder directory, convert to list of full file paths 
    # for ripper jsons
    if type(rippers) != list: 

        rippers = [
            os.path.join(rippers, ripper) for ripper in os.listdir(ripper)
            if ripper.endswith('.json')
        ]

    # initate list to store paths to directories of output data for each 
    # data set 
    write_folders = []

    # iterate through ripper .json files
    for ripper in rippers:

        # create subfolder for this data set 
        write_folder = os.path.join(
            output_folder, 
            os.path.basename(ripper).replace('.json', '')
        )

        print(f'write folder = {write_folder}')

        write_folder = os.path.join(write_folder, 'extracted')

        # make output folder for extracted data if not there already 
        if not os.path.exists(write_folder):
            os.makedirs(write_folder)
            
        # add path to directory of write folder to write_folders
        write_folders.append(write_folder)

        if not os.path.exists(write_folder):
            os.mkdir(write_folder)

        # open ripper .json file as dict
        ripper_dict = open_json(ripper)

        print(f'getting MS1 EICs for {len(compositional_silico_dict)} compositions')
       
        # get MS1 EICs for all compositions
        MS1_EICs = extract_MS1(
            silico_dict=compositional_silico_dict,
            extractor_parameters=extractor_parameters,
            ripper_dict=ripper_dict
        )

        # save MS1 EICs
        write_MS1_EIC_file(
            input_data_file=ripper,
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
            silico_parameters=silico_parameters,
            uniques=False
        )
       
        # write full in silico dict to .json file
        write_pre_fragment_screen_sequence_JSON(
            input_data_file=ripper,
            output_folder=write_folder,
            MSMS_silico_dict=full_MSMS_silico_dict
        )

        print(f'confirming fragments for {len(full_MSMS_silico_dict)} sequences')

        # confirm fragments from silico dict and ripper_dict
        confirmed_fragment_dict = confirm_fragments_sequence_dict(
            silico_dict=full_MSMS_silico_dict,
            ripper_dict=ripper_dict,
            extractor_parameters=extractor_parameters
        )

        print(f'confirmed fragment dict = {confirmed_fragment_dict}')

        # write confirmed fragment silico dict to .json file
        write_confirmed_fragment_dict(
            input_data_file=ripper,
            output_folder=write_folder,
            confirmed_fragment_dict=confirmed_fragment_dict
        )

    # return list of paths to extracted data directories 
    return write_folders 

def standard_postprocess(
    extracted_data_folder,
    postprocess_parameters
):
    
    print(f'extracted_data_folder = {extracted_data_folder}')

    # retrieve full silico MSMS dict from extracted data
    pre_screening_silico_dict = open_json(
        filepath=[
            os.path.join(extracted_data_folder, d_file)
            for d_file in os.listdir(extracted_data_folder)
            if d_file.find('PRE_fragment_screening') > -1][0]
    )

    # retrieve confirmed fragment dict from extracted data 
    confirmed_fragment_dict = open_json(
        filepath=[
            os.path.join(extracted_data_folder, d_file)
            for d_file in os.listdir(extracted_data_folder)
            if d_file.find('confirmed_fragment_dict.json') > -1][0]
    )

    ssw = postprocess_parameters['subsequence_weight']

    if type(ssw) != list: 

        ssw = [ssw]
    
    for subseq_weight in ssw:

        print(f'postprocessing for {subseq_weight} ssw')

        postprocess_params = {
            key: postprocess_parameters[key]
            for key in postprocess_parameters
        }
        postprocess_params['subsequence_weight'] = subseq_weight

        # workout confidence for confirmed sequences at subsequence weight
        confidence_scores = assign_confidence_sequences(
            silico_dict=pre_screening_silico_dict,
            confirmed_dict=confirmed_fragment_dict,
            postprocess_params=postprocess_params
        )
        
        # create new folder directory for current subsequence weight
        confidence_dir = os.path.join(
            os.path.dirname(extracted_data_folder), 
            f'confidence_assignments_{subseq_weight}ssw'
        )

        if not os.path.exists(confidence_dir):
            os.makedirs(confidence_dir)

        write_standard_postprocess_data(
            output_folder=confidence_dir,
            silico_dict=pre_screening_silico_dict,
            confidence_scores=confidence_scores,
            confidence_limit=postprocess_parameters['min_viable_confidence'],
            subsequence_weight=subseq_weight,
            confirmed_fragdict=confirmed_fragment_dict
        )

launch_screen(sys.argv[1])
