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

    # generate compositional silico dict for screening MS1 EICs
    compositional_silico_dict = generate_MS1_compositional_dict(silico_params)

    # if data extraction set to False, make sure not to accidentally call 
    # data filters 
    if not parameters_dict['data_extraction']:
        extractor_params['pre_run_filter'] = False 

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
            if not os.path.isfile(subdir) 
        ]

        extracted_data_dirs = [
            os.path.join(data_dir, 'extracted')
            for data_dir in extracted_data_dirs
        ]
        
        extracted_data_dirs = list(
            filter(
                lambda x: x.find('run_parameters') == -1, 
                extracted_data_dirs
            )
        )

        print(f'extracted data dirs = {extracted_data_dirs}')
    
    if parameters_dict['postprocess']:

        for extracted_data in extracted_data_dirs:

            standard_postprocess(
                extracted_data_folder=extracted_data,
                postprocess_parameters=postprocess_parameters
            )

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
        output_json=os.path.join(output_folder, 'run_parameters.json')
    )

    # apply any pre-filtering steps to ripper data 
    filtered_rippers = apply_filters(
        extractor_parameters=extractor_params,
        ripper_folder=ripper_folder
    )

    # generate dictionary of MS2 signature ion masses for each monomer
    monomer_list = silico_params["MS1"]["monomers"]

    for ripper_json in filtered_rippers:

        # find ripper name
        start = ripper_json.find('filtered_rippers') + 17
        end = ripper_json.find('.json', start)
        ripper_name = ripper_json[start:end]

        print(f'searching {ripper_name} for {monomer_list}')
        ripper_dict = open_json(ripper_json)

        # extract base peak chromatogram (bpc) from each ms1 spectrum
        bpc_dict = {}

        # initiate dictionary for all confirmed monomers found in ripper json
        confirmed_monomers_dict = {}

        # make list of highest peak masses
        highest_peak_masses = []

        for ms1_spectrum_id, ms1_spectrum in ripper_dict["ms1"].items():
            if ms1_spectrum_id.startswith("spectrum_"):

                # find highest peak and it's mass in each spectrum
                # save to bpc dict
                highest_peak = max([
                    float(intensity) for intensity in ms1_spectrum.values() 
                    if type(intensity) != list])

                highest_peak_mass = [
                    mass for mass, intensity in ms1_spectrum.items() 
                    if intensity==highest_peak]

                bpc_dict[ms1_spectrum_id] = [highest_peak_mass[0], highest_peak]

                if highest_peak_mass not in highest_peak_masses:
                    highest_peak_masses.append(highest_peak_mass[0])
        
        # find ms2 spectra that precursors match the highest peak
        for ms2_spectrum_id, ms2_spectrum in ripper_dict["ms2"].items():

            if ms2_spectrum["parent"] in highest_peak_masses:

                # list signature ion dict types
                signature_ion_types = list(MS2_SIGNATURE_IONS.keys())
                signature_ion_type = signature_ion_types[0]

                # in the ms2 spectrum look for the signature ion of the specified monomers
                confirmed_monomers = find_ms2_signature_ions(
                    monomer_list=monomer_list,
                    signature_ion_type=signature_ion_type,
                    spectrum=ms2_spectrum,
                    error=extractor_params["error"],
                    error_abs=extractor_params["err_abs"],
                    signature_ion_dict=MS2_SIGNATURE_IONS)
                        
                if confirmed_monomers:
                    print(f'MS2 signature ions found for {confirmed_monomers} in {ms2_spectrum_id}')

                unconfirmed_monomers = [
                    monomer for monomer in monomer_list
                    if (monomer not in confirmed_monomers
                    and monomer not in MS2_SIGNATURE_IONS["dominant"])
                ]

                unconfirmed_monomers = list(set(unconfirmed_monomers))
                
                # search for unconfirmed monomers
                monomer_mass_dict = {}

                for unconfirmed_monomer in unconfirmed_monomers:

                    # get monomer masses and loss product monomer masses
                    monomer_mass = MONOMERS[unconfirmed_monomer][0]
                                        
                    if unconfirmed_monomer in LOSS_PRODUCTS:
                        loss_products = LOSS_PRODUCTS[unconfirmed_monomer]
                        monomer_masses = [monomer_mass-MASS_DIFF]

                        # get masses - loss products to search for
                        for loss_product in loss_products:
                            monomer_loss_mass = monomer_mass-FG[loss_product]-MASS_DIFF
                            monomer_masses.append(monomer_loss_mass)

                            # add all masses to monomer mass dict
                            monomer_mass_dict[unconfirmed_monomer] = monomer_masses
                                            
                    else:
                        monomer_mass_dict[unconfirmed_monomer] = [(monomer_mass-MASS_DIFF)]

                # search for unconfirmed monomers by looking for mass differences
                # between MS2 peaks
                mass_difference_search = filter_mass_difference(
                    msdata=ms2_spectrum,
                    monomer_massdiffs=monomer_mass_dict,
                    total_comparisons=50,
                    err=extractor_params["error"],
                    err_abs=extractor_params["err_abs"],
                    single_spectrum=True,
                    single_spectra_id=ms2_spectrum_id,
                    ms_level=2)

                # update confirmed monomers dict with mass difference search findings
                monomers_found = []
                for found_monomer_list in mass_difference_search.values():
                    for found_monomer in found_monomer_list:
                        if found_monomer not in monomers_found:
                            monomers_found.append(found_monomer)
                if monomers_found:
                    print(f'{monomers_found} in {ms2_spectrum_id} during mass difference screening')
                # add all confirmed monomers (from ms2 signature screening and
                # mass difference screening) to final confirmed monomers dict
                if confirmed_monomers or monomers_found:
                    confirmed_monomers_dict[ms2_spectrum_id] = {
                        "precursor_mass": ms2_spectrum["parent"],
                        signature_ion_type: confirmed_monomers,
                        "mass_diff": monomers_found
                        }
                    
        # save json for each ripper containing dictionary of format:
        # {"spectrum_id": {"signature_type": [monomers found], "mass_diff": [monomers found]}
        write_to_json(confirmed_monomers_dict, os.path.join(output_folder, f'{ripper_name}_mass_difference_screen.json'))

    return print(f'jsons saved to {output_folder}')


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

        MS1_pre_screening_file = os.path.join(write_folder, "MS1_pre_Screening.json")
    
        write_to_json(
            write_dict=compositional_silico_dict,
            output_json=MS1_pre_screening_file
        )

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
