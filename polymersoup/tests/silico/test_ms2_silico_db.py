'''
This file contains tests for silico module helper functions. These are basic
functions required to perform simple tasks in the silico module, such as
calculate sequence masses or manipulate sequence strings.
'''

import os
import sys
import copy
import json
import pytest
from pymongo import MongoClient

from ...utils.parameter_handlers.parameter_handlers import generate_parameters
from ...silico.polymer_info.polymer import Polymer
from ...silico.helpers.helpers import generate_all_sequences
from ...silico.silico_handler import get_ms2_silico_dict_from_compositions_db

HERE = os.path.abspath(os.path.dirname(__file__))
INPUTS_FOLDER = os.path.join(HERE, '..', 'input_param_files')
SILICO_FOLDER = os.path.join(HERE, '..', 'silico', 'test_silico_dicts')

if sys.platform.startswith('linux'):
    mongo_str = 'mongo'
else:
    mongo_str = 'localhost'

connection = MongoClient(mongo_str, 27017)
connection.drop_database('polymersoup')
polymersoupdb = connection["polymersoup"]


@pytest.fixture
def params():
    '''
    Parameters object for standard peptides.

    Returns:
        Parameters: Parameters object
    '''
    return generate_parameters(
        params_json=os.path.join(
            INPUTS_FOLDER,
            'full_input_parameters_large_sequences.json'))

@pytest.fixture
def depsi_params():
    '''
    Parameters object for depsipeptides.

    Returns:
        Parameters: Parameters object.
    '''
    return generate_parameters(
        params_json=os.path.join(
            INPUTS_FOLDER,
            'full_input_parameters_no_defaults_FGg.json'))

@pytest.fixture
def polymer(params):
    return Polymer(params_obj=params)

@pytest.fixture
def depsi_polymer(depsi_params):
    return Polymer(params_obj=depsi_params)

@pytest.fixture
def polyester_polymer(polyester_params):
    return Polymer(params_obj=polyester_params)

@pytest.fixture
def test_depsi_silico_dict():
    depsi_file = os.path.join(
        SILICO_FOLDER,
        'depsi_test_FGg.json')
    with open(depsi_file) as f:
        load = json.load(f)
        return load

@pytest.fixture
def test_peptide_silico_dict():
    peptide_file = os.path.join(
        SILICO_FOLDER,
        'standard_seqs.json'
    )
    with open(peptide_file) as f:
        load = json.load(f)
        return load

@pytest.fixture
def test_sequences(params, polymer):
    return generate_all_sequences(
        polymer=polymer,
        params=params,
        sequencing=True)

@pytest.fixture
def test_full_ms2_depsi_dict():
    fp = os.path.join(SILICO_FOLDER, 'full_peptide_MSMS_silico.json')
    with open(fp) as f:
        load = json.load(f)
        return load

@pytest.mark.unit
def test_full_fragment_dict_depsi(
    depsi_params,
    depsi_polymer,
    test_full_ms2_depsi_dict
):
    '''
    Tests full MS2 fragment builder for linear depsipeptide sequences: a, b, y
    fragments + signature and composition tags.

    Args:
        depsi_params (Parameters): parameters object (fixture).
        depsi_polymer (Polymer): polymer object (fixture).
        test_full_ms2_depsi_dict (Dict[str, dict]): dict of depsipeptide
            sequences and full MS2 fragments (fixture).
    '''

    local_params = copy.deepcopy(depsi_params)
    local_params.silico.isomeric_targets = [
        seq for seq in test_full_ms2_depsi_dict]

    compositions = list(
        set([''.join(sorted(s)) for s in local_params.silico.isomeric_targets]))

    get_ms2_silico_dict_from_compositions_db(
        params=local_params,
        polymer=depsi_polymer,
        ms1_hits=compositions,
        connection=connection,
        ripper_name='depsi_fragment_test',
        out_folder=SILICO_FOLDER
    )

    for entry in polymersoupdb['depsi_fragment_test_ms2_silico'].find():
        seq = entry['_id']
        frags = entry['series']
        for frag, info in frags.items():
            assert info == test_full_ms2_depsi_dict[seq][frag]
