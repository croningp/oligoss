import os
import json
import logging

# configure info logging
logging.basicConfig(
    filename="logger.log",
    format='%(message)s - %(asctime)s',
    datefmt='%H:%M:%S %m/%d/%Y ',
    level=logging.INFO)

def open_json(filepath):
    """
    Opens parameters json and returns as dictionary.

     Args:
        filepath (str): full file path to parameters json.

    Returns:
        load: parameter dictionary.
    """
    with open(filepath) as f:
        load = json.load(f)
        return load

def write_to_json(
    write_dict,
    output_json
):
    """
    Writes a dict to a json.

    Args:
        write_dict (dict): dictionary to be dumped in json.
        filepath (str): full file path to output json.
    """
    with open(output_json, 'w') as fp:
        json.dump(write_dict, fp, indent=4)

    return logging.info(f'data written to {output_json}')

def return_jsons(input_folder):
    """
    This function creates a list of all json files in the input folder.

    Arguments:
        input_folder {str} -- filepath to mass spec jsons.
    """
    return [os.path.join(
        input_folder,
        file) for file in os.listdir(input_folder) if file.endswith('.json')]
