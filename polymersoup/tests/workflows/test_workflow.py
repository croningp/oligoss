"""
This file contains tests for the full exhaustive screening workflows
"""

import os
import csv
import pytest

from ...utils.parameter_handlers.parameter_handlers import generate_parameters
from ...silico.polymer_info.polymer import Polymer
from ...silico.ms1_silico import generate_ms1_ions
from ...workflows.workflows import (
    exhaustive_screen_forked, exhaustive_screen_spawned)

HERE = os.path.abspath(os.path.dirname(__file__))
INPUTS_FOLDER = os.path.join(HERE, '..', 'input_param_files')
OUTPUTS_FOLDER = os.path.join(HERE, '..', 'output')
RIPPER_FOLDER = os.path.join(HERE, '..', 'extractors', 'test_rippers')

@pytest.fixture
def test_params():
    return generate_parameters(
        os.path.join(INPUTS_FOLDER, 'workflow_test_input_parameters.json'))

@pytest.fixture
def test_polymer(test_params):
    return Polymer(test_params)

@pytest.mark.unit
def compare_forked_spawned_workflows(test_params, test_polymer):
    """
    Compares spawned and forked multiprocessed exhaustive screens to ensure the
    results are identical.

    Args:
        test_params (Parameters): Parameters object.
        test_polymer (Polymer): Polymer object.
    """

    #  generate ms1 ions with three compositions
    ms1_dict = generate_ms1_ions(
        params=test_params,
        polymer=test_polymer,
        sequencing=False
    )

    #  retrieve ripper file
    ripper_file = os.path.join(
        RIPPER_FOLDER,
        '3_compositions_GGGGV_FGGGG_AFGGV.json')

    #  set different outputs for forked and spawned workflows
    forked_output = os.path.join(RIPPER_FOLDER, "forked")
    spawned_output = os.path.join(RIPPER_FOLDER, "spawned")

    #  perform forked screen
    exhaustive_screen_forked(
        ripper=ripper_file,
        compositional_ms1_dict=ms1_dict,
        polymer=test_polymer,
        params=test_params,
        ripper_output=forked_output
    )

    #  perform spawned screen
    exhaustive_screen_spawned(
        ripper=ripper_file,
        compositional_ms1_dict=ms1_dict,
        params=test_params,
        polymer=test_polymer,
        ripper_output=spawned_output
    )

    #  open summary data for forked csv
    forked_csv = [os.path.abspath(f) for f in os.listdir(forked_output)][0]

    #  open summary data for spawned csv
    spawned_csv = [os.path.abspath(f) for f in os.listdir(spawned_output)][0]

    #  make sure results of forked and spawned screen are identical
    assert compare_summary_csv(csvs=[forked_csv, spawned_csv])

def compare_summary_csv(csvs):
    """
    Compares postprocessed summary csvs to make sure sequence assignments are
    identical.

    Args:
        csvs (List[str]): list of full filepaths to target csv files.

    Raises:
        Exception: raised if sequence is missing from one or more csv file.
        Exception: raised if one or more sequence has inconsistent results.

    Returns:
        True: successful run returns True.
    """
    #  init dict to store results for postprocessing
    seq_dict = {}

    for i, f in enumerate(csvs):
        with open(f, "r") as o:
            reader = csv.reader(o, delimiter=",")
            for j, row in enumerate(reader):
                if j > 0:
                    if i == 0:
                        seq_dict[row[0]] == row[1::]
                    else:
                        if row[0] not in seq_dict:
                            raise Exception(
                                f"{row[0]} not confirmed in all runs")
                        if seq_dict[row[0]] != row[1::]:
                            raise Exception(
                                f"data not consistent for {row[0]}")
    return True
