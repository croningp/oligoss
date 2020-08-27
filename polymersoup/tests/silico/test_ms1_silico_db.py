"""
This file contains tests for functions in the ms1_silico module with MongoDB.
"""

import os
import sys
import pytest
from pymongo import MongoClient

from ...utils.parameter_handlers.parameter_handlers import generate_parameters
from ...silico.polymer_info.polymer import Polymer
from ...silico.ms1_silico import generate_ms1_ions_db
from ...silico.helpers.helpers import (
    generate_all_sequences,
    ionize_sequence_precursors
)

HERE = os.path.abspath(os.path.dirname(__file__))
INPUTS_FOLDER = os.path.join(HERE, '..', 'input_param_files')

if sys.platform.startswith('linux'):
    mongo_str = 'mongo'
else:
    mongo_str = 'localhost'
connection = MongoClient(mongo_str, 27017)
connection.drop_database('polymersoup')

@pytest.fixture
def params():
    return generate_parameters(
        params_json=os.path.join(
            INPUTS_FOLDER,
            'full_input_parameters_no_defaults.json'))

@pytest.fixture
def polymer(params):
    return Polymer(params_obj=params)

@pytest.mark.unit
def test_generate_ms1_ions_db(params, polymer):
    """
    Tests MS1 precursor dicts.

    Args:
        params (Parameters): Parameters object
        polymer (Polymer): Polymer object
    """

    Test_1 = False

    silicodb = connection.polymersoup.ms1_silico

    #  generate list of sequences
    sequences = generate_all_sequences(
        polymer=polymer,
        params=params,
        sequencing=True)

    #  generate dict of sequences and corresponding MS1 precursor ions
    generate_ms1_ions_db(
        params=params,
        polymer=polymer,
        sequencing=True,
        connection=connection)

    if silicodb.count() == len(sequences):
        Test_1 = True

    #  check that sequence dict has all target sequences
    assert Test_1 is True

    #  iterate through sequences, making sure precursor m/z values are correct
    for sequence in sequences:

        precursors = silicodb.find_one(
            {'_id': sequence})['precursors']

        calculated_precursors = ionize_sequence_precursors(
            sequence=sequence,
            params=params,
            polymer=polymer)

        assert precursors == calculated_precursors
