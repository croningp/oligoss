"""
This file contains very simple functions for opening, reading from and writing
to files.

"""

import os
import csv
import json
import logging
import multiprocessing

from bson.json_util import dumps, default

from ..utils.run_utils import exception_handler

manager = multiprocessing.Manager()

@exception_handler(fatal=True, verbose=True)
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

@exception_handler(verbose=True)
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

@exception_handler()
def return_jsons(input_folder):
    """
    This function creates a list of all json files in the input folder.

    Arguments:
        input_folder {str} -- filepath to mass spec jsons.
    """
    return [os.path.join(
        input_folder,
        file) for file in os.listdir(input_folder) if file.endswith('.json')]

@exception_handler(fatal=True)
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

@exception_handler(fatal=True)
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

@exception_handler(fatal=True)
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

@exception_handler(timed=True)
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

def return_locks(lock_strs=["csv", "json", "summary"]):
    """
    Get locks for accessing shared outputs.

    Args:
        lock_strs (List[str], optional): string ids for specific locks.
            Defaults to ["csv", "json", "summary"].

    Returns:
        Dict[str, Lock]: dict of string ids and corresponding Lock objects
    """
    return {lock_str: multiprocessing.Lock() for lock_str in lock_strs}

class FileQueueManager():
    def __init__(
        self,
        output_dir,
        sequence_space,
        shared_space_limit=2E7,
        headers=[
            "sequence",
            "confidence",
            "confirmed_core",
            "confirmed_signatures",
            "max_intensity",
            "composition"
        ],
        write_limit=1E5,
        max_queue_size=500
    ):
        """
        Keeps track of queues for writing output data to csv, JSONs in
        multiprocessed workflows.

        Args:
            output_dir (str): Output directory for saving sequencing data.
            sequence_space (int): Total number of sequences to be searched in
                screening full sequence space of products.
            shared_space_limit (Optional[int]): limit on sequence space to write
                results to shared file. If sequence space to screen is above
                this limit, results for each composition will be written to
                their own file. If sequence space is below this limit, results
                will be written to shared output files. Defaults to 2E7.
            headers (Optional[List]): Headers to write on first line of output
                csv files.
            write_limit (Optional[int]): Maximum number of sequences to write to
                single csv file. Defaults to 8E4.
            update_chunk_size (Optional[int]): Maximum number of items to add
                to queue in a single self.csv_queue.put() operation.
            maximum_queue_size(Optional[int]): maximum size of queue before it
                is emptied and contents are written to output csv.

        Properties:
            total_write_count (int): total number of sequences written to output
                sequencing summary files.
            csv_count (int): current iteration of summary csvs.
            write_limit (int): maximum number of sequences to write per file.
            output_dir (str): full path to output directory for dumping data.
            headers (List[str]): headers for output summary csvs.
            csv_queue (Queue): Queue object for controlling write access to
                file.
            output_csv (str): full filepath to current output csv for writing
                data.
            max_q_size (int): maximum size of queue before enforcing write to
                output csv.
            composition_write (bool): specifies whether to write data to
                composition-specific csv or shared output csv file.
        """

        #  keep track of total sequences written to file and current csv number
        #  each csv will contain a maximum of write_limit sequences.
        self.total_write_count = 0
        self.csv_count = 1
        self.write_limit = write_limit
        self.max_q_size = max_queue_size

        if sequence_space >= shared_space_limit:
            self.composition_write = True
        else:
            self.composition_write = False

        #  set paths for output files to write data
        self.output_dir = output_dir
        self.headers = headers
        self.write_limit = write_limit
        self.csv_queue = manager.Queue()
        if not self.composition_write:
            self.write_new_csv()

    def __str__(self):
        msg = (
            f"{self.total_write_count} written to {self.csv_count} CSV files "
            f"so far. Write limit per CSV = {self.write_limit}. Maximum queue"
            f"size set to {self.max_q_size}.")
        if self.composition_write:
            msg += (
                "Data for each composition will be written to a separate file.")
        elif self.csv_queue.empty():
            msg += "Queue is currently empty."
        else:
            msg += f"Queue has approximately {len(self)} items waiting."

        msg += f"{self.total_write_count} sequences written to file so far."

        return msg

    def __len__(self):
        """
        NOTE: the value returned for this is approximate due to delays in
        updating between processes. Use with caution!!!!

        Returns:
            int: number of items remaining in csv_queue (approximate).
        """
        return self.csv_queue.qsize()

    @property
    def output_csv(self):
        return os.path.join(
            self.output_dir, f"sequencing_summary_{self.csv_count}.csv")

    @exception_handler(verbose=True)
    def write_new_csv(self, destination=None):
        """
        Write new output csv.
        """
        if not destination:
            destination = self.output_csv

        if not os.path.exists(self.output_csv):
            with open(self.output_csv, "w") as w:
                writer = csv.writer(w)
                writer.writerow(self.headers)
            logging.info(f"output csv written to {self.output_csv}")

    def queue_data(self, data):
        if not self.composition_write:
            self.csv_queue.put(data)
        else:
            composition = data[-1]
            composition_csv = os.path.join(
                self.output_dir, f"sequencing_summary_{composition}.csv")
            self.write_new_csv(destination=composition_csv)
            with open(composition_csv, "a") as a:
                self.total_write_count += 1
                writer = csv.writer(a)
                writer.writerow(data)

    def update(self):
        """
        Updates csv_queue and writes sequencing data to csv.

        """
        #  as long as there is stuff in the queue, keep writing it to file
        while not self.csv_queue.empty():

            #  get latest item in queue
            data = self.csv_queue.get()

            #  ensure current output_csv hasn't exceeded write limit
            if self.total_write_count < self.csv_count * self.write_limit:
                self.append_csv(data)
                self.total_write_count += 1
            else:

                #  write limit has been exceeded for current output_csv
                #  => move on to new output_csv, updating csv_count and
                #  adding remaining unwritten data back to queue
                self.csv_count += 1
                self.write_new_csv()
                self.csv_queue.put(data)

            self.csv_queue.task_done()

    @exception_handler(verbose=False)
    def append_csv(self, data):
        with open(self.output_csv, "a") as w:
            writer = csv.writer(w)
            writer.writerow(data)
