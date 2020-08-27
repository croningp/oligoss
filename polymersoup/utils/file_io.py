"""
This file contains very simple functions for opening, reading from and writing
to files.

"""

import os
import csv
import json
from bson.json_util import dumps, default

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
    output_json,
    type=dict
):
    """
    Writes a dict to a json.

    Args:
        write_dict (dict): dictionary to be dumped in json.
        filepath (str): full file path to output json.
    """
    with open(output_json, 'w') as fp:
        if type == dict:
            json.dump(write_dict, fp, indent=4)
        elif type == str:
            json.dumps(obj=write_dict)

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

def append_locked_json(fpath, lock, dump_dict):
    """
    Appends contents of a dict to a new line in JSON as string.

    Args:
        fpath (str): full file path to output JSON.
        lock (Lock): Lock object.
        dump_dict (dict): dict to dump as string.
    """
    #  acquire lock, open and dump dict to JSON, release lock
    lock.acquire()

    #   open JSON file for reading and writing
    with open(fpath, mode="r+") as o:

        #  go to end of file, set position at end of file so dict can be
        #  dicted in new line
        o.seek(0, 2)
        position = o.tell()
        o.seek(position)

        #  dump dict as string on a new line
        o.write(f"{json.dumps(dump_dict)}")

    lock.release()

def mongodb_to_json(
    collection_name,
    output_json
):
    """ This function dumps a single mongoDB collection to a json file.

    Args:
        collection_name (str): database collection name.
        output_json (str): full filepath of output json.
    """

    cursor = collection_name.find({})

    with open(output_json, 'w') as file:
        file.write('[')
        doc_count = 0
        for document in cursor:

            # keep count of documents
            doc_count += 1
            num_max = cursor.count()

            if (num_max >= 1) and (doc_count <= num_max - 1):
                dump = dumps(
                    document, sort_keys=False, indent=4, default=default)
                file.write(dump)
                file.write(',')

            # if document is last, don't add comma
            elif (num_max == doc_count):
                dump = dumps(
                    document, sort_keys=False, indent=4, default=default)
                file.write(dump)

        file.write(']')
