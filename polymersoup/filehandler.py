"""

This file contains functions for reading and passing on run parameters
for de novo peptide polymer soup sequencing

"""
import os
import csv
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
        json.dump(write_dict, fp, indent=4)

    print(f'data written to  {output_json}')

def write_new_csv(
    csv_file,
    headers=None,
    delimitingchar=','
):
    """
    Writes a new .csv file

    Args:
        csv_file (str): full file path to new .csv file to be written
        headers (list, optional): list of headers to be written to first
            row of .csv file. Defaults to None.
        delimitingchar (str, optional): delimiter for separating columns in
            output .csv file. Defaults to ','.
    """
    with open (csv_file, 'w') as write_file:

        write = csv.writer(write_file, delimiter = delimitingchar)
        if headers:
            write.writerow(headers)

def append_csv_row(
    csv_file,
    append_list,
    delimitchar=','
):
    """
    Appends a new row to existing .csv file

    Args:
        csv_file (str): full file path to pre-written .csv file
        append_list (list): list of strings to write to new .csv row
        delimitchar (str, optional): delimiter for separating columns in
            .csv file. Defaults to ','.
    """
    with open(csv_file, 'a') as ofile:
        writer = csv.writer(ofile, delimiter=delimitchar)
        writer.writerow(append_list)


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

def write_EIC_file(
    input_data_file,
    output_folder,
    EICs,
    ms_level=1
):
    """
    Writes an EIC dict to a .json file

    Args:
        input_data_file (str): full file path to mzml ripper data file used to
            generate MS1 EICs
        output_folder (str): path to output folder directory, where data will
            be saved
        MS1_EICs (dict): dictionary of sequences and / or compositions and
            their corresponding MS1 EICs
        ms_level (int, optional): specifies MS level of EIC data. Defaults to
            1
    """

    input_file = os.path.basename(input_data_file).replace(".json", "")

    output_file = os.path.join(output_folder, f'{input_file}_MS{ms_level}_EICs.json')

    write_to_json(
        write_dict=EICs,
        output_json=output_file
    )
    print(f'MS{ms_level} EICs written to {output_file}')

def write_pre_fragment_screen_sequence_JSON(
    input_data_file,
    output_folder,
    MSMS_silico_dict
):
    """
    Writes in silico sequence dict for sequences that have been detected by
    MS1 composition - i.e. BEFORE CONFIRMING MS2 FRAGMENTS AND ASSIGNING
    CONFIDENCE

    Args:
        input_data_file (str): full file path to input mzml ripper data file
        output_folder (str): directory of output folder, where data will be
            saved
        MSMS_silico_dict (dict): in silico sequence dict for sequences,
            potential MS1 precursors and MS2 fragments
    """
    input_file = os.path.basename(input_data_file).replace(".json", "")
    input_file = input_file + "_PRE_fragment_screening_silico_dict.json"

    output_file = os.path.join(output_folder, input_file)

    write_to_json(
        write_dict=MSMS_silico_dict,
        output_json=output_file
    )

    print(f'full MSMS silico dict for sequences detected at MS1 written to:')
    print(output_file)

def write_confidence_assignments(
    input_data_file,
    output_folder,
    confidence_assignments
):
    """
    Writes confidence assignments for sequences with confirmed MS1 precursors
    and MS2 fragments to a .json file

    Args:
        input_data_file (str): full file path to input mzml ripper data file
        output_folder (str): directory of output folder, where data will be
            saved
        confidence_assignments (dict): dictionary of sequences and their
            assigned confidence scores
    """

    input_file = os.path.basename(input_data_file).replace(".json", "")
    input_file = input_file + "_sequence_confidence_assignments.json"

    output_file = os.path.join(output_folder, input_file)

    write_to_json(
        write_dict=confidence_assignments,
        output_json=output_file
    )

    print(f'confidence assignments written to {output_file}')

def write_confirmed_fragment_dict(
    input_data_file,
    output_folder,
    confirmed_fragment_dict
):
    """
    [summary]

    Args:
        input_data_file (str): full file path to input mzml ripper data file
        output_folder (str): directory of output folder, where data will be
            saved
        confirmed_fragment_dict (dict): dictionary of sequences and confirmed
            fragments for all sequences that have ANY confirmed MS2 fragments
    """
    input_file = os.path.basename(input_data_file).replace(".json", "")
    input_file = input_file + "_confirmed_fragment_dict.json"

    output_file = os.path.join(output_folder, input_file)

    write_to_json(
        write_dict=confirmed_fragment_dict,
        output_json=output_file
    )

    print(f'confirmed fragment dict written to {confirmed_fragment_dict}')

def write_unique_fragment_dict(
    input_data_file,
    output_folder,
    unique_fragment_dict
):
    """
    Writes a dictionary of confidently assigned sequences and their
    corresponding in silico fragment data for fragments that have been confirmed,
    including unique fragments

    Args:
        input_data_file (str): full file path to input mzml ripper data file
        output_folder (str): directory of output folder, where data will be
            saved
        unique_fragment_dict (dict): dictionary of confident sequences and
            corresponding in silico data for confirmed fragments, including
            unique confirmed fragments
    """
    input_file = os.path.basename(input_data_file).replace(".json", "")
    input_file = input_file + "_confirmed_fragments_with_uniques.json"

    output_file = os.path.join(output_folder, input_file)

    write_to_json(
        write_dict=unique_fragment_dict,
        output_json=output_file
    )

    print(f"unique fragment dictionary written to {output_file}")

def write_standard_postprocess_data(
    output_folder,
    silico_dict,
    confirmed_fragdict, 
    confidence_scores,
    confidence_limit,
    subsequence_weight,
    MS1_EICs
):
   
    write_to_json(
        write_dict=confidence_scores,
        output_json=os.path.join(output_folder, 'confidence_scores.json')
    )

    
    confirmed_fragdict = {
        seq: {
            "MS1": silico_dict[seq]["MS1"],
            "MS2": {
                frag: masses 
                for frag, masses in silico_dict[seq]["MS2"].items()
                if (frag in confirmed_fragdict[seq]
                and frag != "signatures")
            },
            "confirmed_signatures": [
                frag for frag in confirmed_fragdict[seq]
                if frag in silico_dict[seq]["MS2"]["signatures"]
            ],
            "peak_list": silico_dict[seq]["peak_list"]
        }
            for seq in confirmed_fragdict
    }

    print(f'confirmed_fragdict = {confirmed_fragdict}')

    write_to_json(
        write_dict=confirmed_fragdict,
        output_json=os.path.join(
            output_folder,
            f"confirmed_frag_dict.json"
        )
    )

    summary_csv = os.path.join(
        output_folder, 
        f"postprocess_summary_{subsequence_weight}ssw.csv"

    )
    write_new_csv(
        csv_file=summary_csv,
        headers=[
            "Sequence", 
            "Confidence_Score", 
            "Confirmed_Fragments", 
            "Confirmed_signatures",
            "Max_Intensity"
        ]
    )
    
    # get max peak intensity of each sequence from its MS1 EIC
    max_intensities = {
        seq: max([Rt_I[1] for Rt_I in MS1_EICs["".join(sorted(seq))]])
        for seq in confirmed_fragdict
    }

    # write all data to final csv for each sequence 
    for seq in confirmed_fragdict: 

        append_csv_row(
            csv_file=summary_csv,
            append_list=[
                seq,
                confidence_scores[seq],
                [frag for frag in confirmed_fragdict[seq]["MS2"]
                if frag != 'signatures'],
                confirmed_fragdict[seq]["confirmed_signatures"],
                max_intensities[seq]
            ]
        )
    

