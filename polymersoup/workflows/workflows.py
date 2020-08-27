"""
This file contains functions for the main Polymersoup experimental workflows
"""
import os
import multiprocessing

#  import silico functions
from ..silico.polymer_info.polymer import Polymer

#  import extractors functions
from ..extractors.filters import mzml_to_json, prefilter_all

#  import from utils
from ..utils.file_io import return_jsons
from ..utils.run_utils import RunInfo

#  import from other workflows submodules
from .workflow_helpers.linux_workflows import screen_rippers_linux
from .workflow_helpers.windows_workflows import screen_rippers_windows

#  get lock object for synchronising access to summary csv, silico and spectral
#   assignment JSON files between processes
CSV_LOCK, JSON_LOCK = (
    multiprocessing.Lock(),
    multiprocessing.Lock()
)


def exhaustive_screen_multi(params, data_folder, out_folder):
    """
    Exhaustive screen with multiprocessing

    Args:
        params (Parameters): parameters object.
        data_folder (str): ripper data dir.
        out_folder (str): output data dir
    """

    #  retrieve info for OS and number of cores
    run_info = RunInfo(params=params)
    run_info.log_record()

    #  generate and log Polymer object
    polymer = Polymer(params)
    polymer.log_record()

    #  make sure all out files get dumped in "extracted" subfolder
    out_folder = os.path.join(out_folder, "extracted")
    if not os.path.exists(out_folder):
        os.mkdir(out_folder)

    #  convert mzml files to ripper JSONs
    mzml_to_json(
        input_folder=data_folder,
        extractor_parameters=params.extractors
    )

    #  get list of filepaths to ripper JSON files
    rippers = return_jsons(data_folder)

    #  pre - filter ripper JSONs and save output
    filtered_rippers = prefilter_all(
        rippers=rippers,
        extractor_parameters=params.extractors
    )

    if run_info.init_method == "forked":
        screen_rippers_linux(
            filtered_rippers=filtered_rippers,
            params=params,
            polymer=polymer,
            out_folder=out_folder,
            run_info=run_info)
    else:
        screen_rippers_windows(
            filtered_rippers=filtered_rippers,
            params=params,
            polymer=polymer,
            out_folder=out_folder,
            run_info=run_info
        )
