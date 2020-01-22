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
    """ This function reads an input parameters .json file and decides what
    kind of screening method to use, then performs sequencing screen.
    
    Args:
        input_parameters_file (str) -- full file path to input parameters .json
            file.
    
    Raises:
        Exception: This package requires Python version 3.6 or later.
    """
    
    # check version of python
    if sys.version_info < (3, 6):
        print(f"you are running Python version {sys.version}")
        raise Exception("this package requires Python version 3.6 or later")

    # load parameters
    parameters_dict = generate_parameters_dict(input_parameters_file)

    # perform exhaustive screen if that is specified in input parameters
    if parameters_dict['screening_method'] == 'exhaustive':
        exhaustive_screen(parameters_dict)

    if parameters_dict['screening_method'] == 'mass_difference':
        mass_difference_screen(parameters_dict)

def exhaustive_screen(parameters_dict): 
    """ This function performs an exhaustive or 'brute force' screen for all
    possible sequences arising from input monomers and constraints set, read
    directly from input parameters .json file.
    
    Args:
        parameters_dict (dict) -- input parameters dictionary in standard
            input parameters format.
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
        output_json=os.path.join(output_folder, 'run_parameters.json'))

    # load parameters for postprocessing operations
    postprocess_parameters = parameters_dict["postprocessing_parameters"]

    # load output folder for saving data, create if it does not exist
    output_folder = directories['output_folder']
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)
    
    write_to_json(
        write_dict=parameters_dict,
        output_json=os.path.join(output_folder, 'run_parameters.json'))

    if parameters_dict["perform_silico"]:
        print(f'generating MS1 silico dict')
        # generate compositional silico dict for screening MS1 EICs
        compositional_silico_dict = generate_MS1_compositional_dict(silico_params)
    else:
        
        MS1_s_json = os.path.join(output_folder, 'extracted', 'MS1_pre_Screening.json')
        try: 
            compositional_silico_dict = open_json(MS1_s_json)
        except IOError: 
            raise Exception(f'the file {MS1_s_json} does not exist. Please check whether you already have this data')

    # if data extraction set to False, make sure not to accidentally call 
    # data filters 
    if not parameters_dict['data_extraction']:
        extractor_params['pre_run_filter'] = False 

    # apply any pre-filtering steps to ripper data 
    filtered_rippers = apply_filters(
        extractor_parameters=extractor_params,
        ripper_folder=ripper_folder)

    # if data_extraction set to True, extract data 
    if parameters_dict['data_extraction']:

        extracted_data_dirs = standard_extraction(
            rippers=filtered_rippers,
            output_folder=output_folder,
            silico_parameters=silico_params,
            compositional_silico_dict=compositional_silico_dict,
            extractor_parameters=extractor_params)

    else:

        extracted_data_dirs = [
            os.path.join(output_folder, subdir) 
            for subdir in os.listdir(output_folder)
            if not os.path.isfile(subdir)]

        extracted_data_dirs = [
            os.path.join(data_dir, 'extracted')
            for data_dir in extracted_data_dirs]
        
        extracted_data_dirs = list(
            filter(
                lambda x: x.find('run_parameters') == -1, 
                extracted_data_dirs))

        print(f'extracted data dirs = {extracted_data_dirs}')
    
    if parameters_dict['postprocess']:

        for extracted_data in extracted_data_dirs:

            standard_postprocess(
                extracted_data_folder=extracted_data,
                postprocess_parameters=postprocess_parameters)

def mass_difference_screen(parameters_dict): 
    """ This function screens for mass differences based on input monomers and
    constraints set, read directly from input parameters .json file. Useful for
    identifying monomers present in the data.
    
    Args:
        parameters_dict (dict) -- input parameters dictionary in standard
            input parameters format.

    Returns:
        json (dict) -- {"spectrum_id": {"signature_type": [monomers found],
            "mass_diff": [monomers found]}
    """
    
    # load parameters for in silico operations
    silico_params = parameters_dict['silico_parameters']

    # load parameters for data extraction operations
    extractor_params = parameters_dict['extractor_parameters']

    print(f'{extractor_params.keys()}')
    # load important file paths from parameters_dict
    directories = parameters_dict['directories']

    # load location of mass spec data in mzml_ripper .json file format
    ripper_folder = directories['ripper_folder']

    # load output folder for saving data, create if it does not exist
    output_folder = directories['output_folder']
    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    # write run parameters to .json
    write_to_json(
        write_dict=parameters_dict,
        output_json=os.path.join(output_folder, 'run_parameters.json'))

    # apply any pre-filtering steps to ripper data 
    filtered_rippers = apply_filters(
        extractor_parameters=extractor_params,
        ripper_folder=ripper_folder)

    # generate dictionary of MS2 signature ion masses for each monomer
    monomer_list = silico_params["MS1"]["monomers"]

    # write csv to store summary info for spectral assignments for ALL input
    # samples
    summary_csv = os.path.join(output_folder, f'{monomer_list}_massdiff_screen_summary.csv')
    write_new_csv(
        csv_file=summary_csv,
        headers=[
            'Sample',
            'N_total_spectra', 
            'total_assigned_spectra', 
            '%_assigned_spectra',
            'N_basepeak_assignments',
            '%_basepeak_assignments',
            '%_assigned_precursors',
            'N_assigned_precursors'])

    # iterate through rippers and make monomer-specific spectral fingerprint
    # assignments, generate base peak chromatograms (bpcs)
    for ripper_json in filtered_rippers:
        
        # find ripper name
        start = ripper_json.find('filtered_rippers') + 17
        end = ripper_json.find('.json', start)
        ripper_name = ripper_json[start:end]

        print(f'searching {ripper_name} for {monomer_list}')
        ripper_dict = open_json(ripper_json)

        # generate ms1 bpc to find most intense potential precursors for 
        # screening @ ms2 below 
        ms1_bpc = generate_bpc(
            msdata=ripper_dict,
            N_peaks_per_spectrum=extractor_params["N_bpc_peaks"],
            ms_level=1)
        
        # retrieve monomer massdiff and signature ion fingerprints 
        monomer_fingerprints = generate_monomer_massdiff_signature_fingerprints(
            monomers=monomer_list,
            losses=True,
            signature_types=None)

        amines, aldehydes = ['x', 'p', 'e'], ['h', 'o', 't']

        MS1_dict = generate_ms1_mass_dictionary_adducts_losses(
            monomers=[monomer for monomer in monomer_list if len(monomer)==1],
            adducts=silico_params["MS1"]["ms1_adducts"],
            monomer_modifications_dict={},
            mode='pos',
            min_length=silico_params["MS1"]["min_length"],
            max_length=silico_params["MS1"]["max_length"],
            universal_monomer_modification=False,
            universal_terminal_modification=False,
            min_z=1,
            max_z=None,
            sequencing=False,
            terminal_modifications_dict={},
            max_total_losses=2)
        
        MS1_dict = {
            seq: seq_masses
            for seq, seq_masses in MS1_dict.items()
            if (
                abs(
                    len([c for i, c in enumerate(seq) if c in amines])
                    -len([c for i, c in enumerate(seq) if c in aldehydes]))
                <= 1)
        }

        
        # check whether precursors are to be filtered by bpc
        if not extractor_params["massdiff_bpc_filter"]:
            bpc_dict = None
        else:
            bpc_dict = ms1_bpc
        
        # identify spectra with monomer fingerprints 
        fingerprint_spectra = fingerprint_screen_MSn_from_bpc_precursors(
            bpc_dict=bpc_dict,
            silico_info=monomer_fingerprints,
            err=extractor_params['error'],
            err_abs=extractor_params['err_abs'],
            ms_level=2,
            comparisons_per_spectrum=None,
            msdata=ripper_dict)

        precursor_assignments = screen_MS1_precursors_masslib(
            MS1_masslib=MS1_dict,
            extractor_parameters=extractor_params,
            spectral_assignments=fingerprint_spectra)
        
        if "assigned_precursors" in precursor_assignments:
            N_precursor_assignments = len(precursor_assignments["assigned_precursors"])
            precursor_assignment_pc = (N_precursor_assignments/(len(fingerprint_spectra)-1))*100
        else:
            N_precursor_assignments, precursor_assignment_pc = 0, 0

        write_to_json(
            write_dict=precursor_assignments,
            output_json=os.path.join(
                output_folder,
                f'{ripper_name}_precursor_assignments.json'))

        write_to_json(
            write_dict=MS1_dict,
            output_json=os.path.join(
                output_folder,
                f'{ripper_name}_MS1_lib.json'))

        # write spectral assignments with identified monomer fingerprints to 
        # json file
        write_to_json(
            write_dict=fingerprint_spectra, 
            output_json=os.path.join(
                output_folder, 
                f'{ripper_name}_mass_difference_screen.json'))

        # summarise monomer fingerprint assignments for sample
        summary_info = postprocess_massdiff_spectral_assignments(
            spectral_assignment_dict=fingerprint_spectra,
            ripper_dict=ripper_dict)

        # add summary of fingerprint spectral assignments to summary csv
        append_csv_row(
            csv_file=summary_csv,
            append_list=[
            ripper_name,
            summary_info['N_total_spectra'], 
            summary_info['total_assigned_spectra'], 
            summary_info['%_assigned_spectra'],
            summary_info['N_basepeak_assignments'],
            summary_info['%_basepeak_assignments'],
            precursor_assignment_pc,
            N_precursor_assignments])

        # init string to generat bpc file name and write ms1 bpc to json file
        bpc_string = f'_bpc_{extractor_params["N_bpc_peaks"]}peaks_per_spectrum'
        write_to_json(
            write_dict=ms1_bpc,
            output_json=os.path.join(
                output_folder,
                f'{ripper_name}{bpc_string}.json'))

def apply_filters(
    extractor_parameters,
    ripper_folder
):
    """
    Applies any pre-run filters to ripper data and returns list of full 
    file paths to ripper jsons.
    
    Args:
        extractor_parameters (dict): dictionary of extractor parameters, 
            retrieved from input parameters JSON file.
        ripper_folder (str): full file path to folder containing mzml ripper
            JSON files.
    
    Returns:
        list: list of full file paths to mzml ripper JSON files.
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
    """ This function generates full insilico fragmentation MSMS data for all possible
        sequences that match compositions found in MS1 EICs. It confirms fragments 
        and writes them to a confirmed fragment silico dict json.
    
    Args:
        rippers (str or list): full file path to folder containing mzml ripper
            JSON files, generated by mzml ripper package.
        extractor_parameters (dict): dictionary of extractor parameters, 
            retrieved from input parameters JSON file.
        silico_parameters(): dictionary of silico parameters, 
            retrieved from input parameters JSON file.
        compositonal_silico_dict(dict): MS1 compositional dict.
        output_folder(str): full file path desired location of polymermassspec output.
    
    Returns:
        list -- list of paths to extracted data directories.
    """
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
            os.path.basename(ripper).replace('.json', ''))

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

        MS1_pre_screening_file = os.path.join(write_folder, "MS1_pre_Screening.json")
    
        write_to_json(
            write_dict=compositional_silico_dict,
            output_json=MS1_pre_screening_file)

        # get MS1 EICs for all compositions
        MS1_EICs = extract_MS1(
            silico_dict=compositional_silico_dict,
            extractor_parameters=extractor_parameters,
            ripper_dict=ripper_dict
        )

        # save MS1 EICs
        write_EIC_file(
            input_data_file=ripper,
            output_folder=write_folder,
            EICs=MS1_EICs,
            ms_level=1
        )

        # remove compositions that are not present at sufficient abundance
        # from in silico compositional dict
        compositional_silico_dict = {
            composition: compositional_silico_dict[composition]
            for composition in MS1_EICs
            if composition.find('retention') == -1
        }

        print(f'generating full MSMS dict for {len(MS1_EICs)} compositions')

        # generate full insilico fragmentation MSMS data for all possible
        # sequences that match compositions found in MS1 EICs
        full_MSMS_silico_dict = generate_MSMS_insilico_from_compositions(
            composition_dict=compositional_silico_dict,
            silico_parameters=silico_parameters,
            uniques=False)
       
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
            extractor_parameters=extractor_parameters)

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
    """ This function retrieves extracted MS data and works out confidence
        scores, Rt, I etc. for each confirmed sequence at all subsequence weights.
    
    Arguments:
        extracted_data_folder (str) -- folder where extracted MS data is.
        postprocess_parameters (dict) -- dictionary of postprocessing parameters, 
            retrieved from input parameters JSON file.
    """
    
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

    # retrieve MS1 EICs dict from extracted data 
    MS1_EICs = open_json(
        filepath = [
            os.path.join(extracted_data_folder, d_file)
            for d_file in os.listdir(extracted_data_folder)
            if d_file.find('MS1_EIC') > -1][0]
    )

    # retrieve subsequence weight(s) - the relative weighting for 
    # continuous fragment coverage for confidence assignment
    ssw = postprocess_parameters['subsequence_weight']

    # ensure ssw is handled as a list
    if type(ssw) != list: 
        ssw = [ssw]
    
    # iterate through each subsequence weight and assign a confidence
    # score
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
            confirmed_fragdict=confirmed_fragment_dict,
            MS1_EICs=MS1_EICs
        )

launch_screen(sys.argv[1])
