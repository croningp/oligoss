"""

This file contains functions for reading and passing on run parameters
for de novo peptide polymer soup sequencing

"""
import os
import json

def __init__():
    pass

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

def write_to_json(
    write_dict,
    output_json
):
    """
    Writes a dict to a json
    Args:
        write_dict (dict): dictionary to be dumped in json
        filepath (str): full file path to output json
    """

    with open(output_json, 'w') as fp:
        json.dump(write_dict, fp)

    print(f'data written to  {output_json}')

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

def list_files_by_type(
    folder,
    filetype
):
    """
    Takes a folder and lists full file paths for all files within that folder
    that match a specified type - either file extension or project-specific
    files 'sequence_json' or 'ripper' files

    Args:
        folder (str): full path to folder containing target files
        filetype (str): specifies filetype - either file extension such as
                '.json', '.mzml', '.csv' etc or a project-specific category
                of files 'sequence_json' or 'ripper'

    Returns:
        files (list): list of full file paths for all files in folder that
                match specified filetype
    """
    # checks to see if any project-specific filetypes have been specified -
    # 'ripper' = json files produced from mzml_ripper
    # 'sequence_json' = json files containing in silico sequence libraries
    if filetype == 'ripper':

        # 'ripper' substring is unique to mzml_ripper files
        filetype, substring = '.json', 'ripper'

    elif filetype == 'sequence_json':

        # 'monomers' substring is unique to sequence json files containing
        # in silico sequence libraries
        filetype, substring = '.json', 'monomers'
    else:

        # if no project-specific filetype has been specified, assume filetype
        # is generic file extension with no additional substrings to search for
        filetype, substring = filetype, ""

    files = [
        os.path.join(folder, file)
        for file in os.listdir(folder)
        if file.endswith(filetype)
        and file.find(substring) > -1
    ]

    return files

def generate_insilico_writefile_string(
    folder,
    silico_dict
):
    """
    Creates a full filepath for in silico sequence json

    Args:
        folder (str): path to folder where in silico json is to be written to
        silico_dict (dict): MSMS sequence json from SilicoGenerator

    Returns:
        str: full filepath to in silico sequence json file that is to be
                written and populated with in silico sequence data
    """
    monomers = silico_dict["MS1"]["monomers"]
    mode = silico_dict["mode"]
    ms1_adducts = silico_dict["MS1"]["ms1_adducts"]
    max_length = silico_dict["MS1"]["max_length"]
    min_length = silico_dict["MS1"]["min_length"]
    start_tags = silico_dict["MS1"]["terminal_tags"]["0"]
    end_tags = silico_dict["MS1"]["terminal_tags"]["-1"]

    write_string = f"monomers={monomers},mode={mode},adducts={ms1_adducts},min_len={min_length},max_len={max_length}"
    if start_tags:
        write_string = write_string + f"0terminaltags={start_tags}"
    if end_tags:
        write_string = write_string + f"-1terminaltags={end_tags}"
    write_string = f"{write_string}.json"
    return os.path.join(folder, write_string)


    