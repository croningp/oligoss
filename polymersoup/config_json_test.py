import os
import json

def read_parameters(parameters_file):
    """
    Takes parameter json and returns subsets of parameter info that are
    required for a run

    Args:
        parameters_file (str): full filepath to parameter json, which should
        specify all parameters needed to conduct a full run of in silico
        fragmentation, data extraction, sequence confirmation and
        postprocessing, or any combination of these steps - see Readme
    """
    # open the parameters json as a dictionary to be read
    params_dict = open_json(parameters_file)

    return params_dict

def open_json(filepath):
    """
    Opens parameters json and returns as dictionary

    Args:
        filepath (str): full file path to parameters json

    Returns:
        load: parameter dictionary
    """
    with open(filepath) as f:
        load = json.load(f)
        return load

def return_polymer_config(
    parameters_file
):
    """ This function reads the parameters file for the polymer type
    and returns the relevant polymer config file as a dictionary.
    
    Arguments:
        parameters_file {str} -- filepath for the parameters .json

    Returns:
        dict version of polymer configuration file
    """
    # read parameters file
    parameters_dict = read_parameters(parameters_file)

    # get polymer type from the input parameters .json
    polymer_directories = parameters_dict["directories"]
    
    polymer_type = polymer_directories["polymer_type"]

    # generate polymer config file name
    polymer_file = f'{polymer_type}.json' 
    
    # get full polymermassspec/polymersoup filepath
    polymermassspec_path = os.path.dirname(os.path.abspath(__file__))

    # generate the polymer config file path
    polymer_config_file_path = os.path.join(polymermassspec_path, "insilico", "Config_files", polymer_file)

    # read polymer config file
    polymer_config_dict = open_json(polymer_config_file_path)
    print(polymer_config_dict)
    return polymer_config_dict

print(return_polymer_config("C:/Users/group/polymermassspec/Examples/InputParams_Test.json"))