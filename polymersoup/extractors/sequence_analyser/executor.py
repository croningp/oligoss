"""Run script for mass processing sequencing data

Can convert mzML files in bulk into json data
Parses every MSn data file and searches for fragments in the sequencing data

.. moduleauthor:: Graham Keenan 2019
.. moduleauthor:: David Doran 2019

"""

import os
import json
import filetools as ft
import mzmlripper.extractor as mzml
from sequence_analyser import SequenceAnalyser


def read_json(filename: str) -> dict:
    """Reads a JSON file and returns the data

    Arguments:
        filename {str} -- Name of the JSON file

    Returns:
        dict -- JSON data
    """

    with open(filename) as f_d:
        return json.load(f_d)


def get_json_files(json_dir: str) -> list:
    """Gets all JSON files in a directory

    Arguments:
        json_dir {str} -- Path to the directory

    Returns:
        list -- List of absolute filepaths of all files within
    """

    # Get all JSON files in given directory
    return [
        file for file in ft.list_files(json_dir)
        if file.endswith(".json")
    ]


def process_msdata(
        sequence_filename: str,
        sequence_data: dict,
        ms_files: list,
        output_dir: str
    ):
    """Processes each Mass Spec file to search for sequence information

    Arguments:
        sequence_filename {str} -- Name of the sequence file
        sequence_data {dict} -- Sequencing Information
        ms_files {list} -- List of all Mass Spec JSON files
        output_dir {str} -- Directory to store output
    """

    # Process each Mass Spec file
    for ms in ms_files:
        msdata = read_json(ms)
        ms_name = ms.replace(".json", "")
        filename = sequence_filename + "_" + os.path.basename(ms_name)
        SequenceAnalyser(filename, msdata, sequence_data, output_dir, bruker=True).parse_data()


def run(ms_json_dir: str, seq_json_dir: str, output_dir: str):
    """Runs the sequencing analysis in bulk
    Parses each Sequence file and reads all MS data files to look for data

    Arguments:
        ms_json_dir {str} -- Path to Mass Spec JSON files
        seq_json_dir {str} -- Path to Sequencing JSON files
        output_dir {str} -- Where to place the output of the program

    Keyword Arguments:
        mzml_folder {str} -- Path to mzML files if they need parsed (default: {""})
    """

    # Get all sequence JSON filepaths
    seq_files = get_json_files(seq_json_dir)

    # Get all mass spec JSON filepaths
    msfiles = get_json_files(ms_json_dir)

    # For every sequence file, process the Mass Spec data
    for seq in seq_files:
        seq_filename = seq.split(os.sep)[-1].replace(".json", "")
        seq_data = read_json(seq)
        process_msdata(seq_filename, seq_data, msfiles, output_dir)


if __name__ == "__main__":
    run("ripper", "sequences", "data")
