"""
This file contains very simple functions for opening, reading from and writing
to files.

"""

import os
import csv
import json

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

def return_jsons(input_folder):
    """
    This function creates a list of all json files in the input folder.

    Arguments:
        input_folder {str} -- filepath to mass spec jsons.
    """
    return [os.path.join(
        input_folder,
        file) for file in os.listdir(input_folder) if file.endswith('.json')]

def write_locked_csv(fpath, lock, write_list=None, delimiter=','):
    """
    Writes a new csv.

    Args:
        fpath (str): full filepath to output csv file.
        lock (object): Lock object for synchronising access to csv between
            processes.
        write_list (list, optional): list to write to 1st row of new csv.
            Defaults to None.
        delimiter (str, optional): delimiting char. Defaults to ','.
    """
    lock.acquire()
    with open(fpath, 'w') as o:
        writer = csv.writer(o, delimiter=delimiter)
        if write_list:
            writer.writerow(write_list)
    lock.release()

def append_locked_csv(fpath, lock, write_list, delimiter=','):
    """
    Appends to a new row of a csv.

    Args:
        fpath (str): full filepath to csv.
        lock (object): Lock object for synchronising access to csv between
            processes.
        write_list (list): list to append in new row.
        delimiter (str, optional): delimiting char. Defaults to ','.
    """

    #  acquire lock, open and write to csv, release lock
    lock.acquire()
    with open(fpath, 'a') as o:
        writer = csv.writer(o, delimiter=delimiter)
        writer.writerow(write_list)
    lock.release()
