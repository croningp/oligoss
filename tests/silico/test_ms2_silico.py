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

from oligoss.utils.parameter_handlers.parameter_handlers import generate_parameters
from oligoss.silico.polymer_info.polymer import Polymer
from oligoss.silico.helpers.helpers import generate_all_sequences
from oligoss.silico.ms2_silico import build_linear_fragments_sequence_dict

HERE = os.path.abspath(os.path.dirname(__file__))
INPUTS_FOLDER = os.path.join(HERE, '..', 'input_param_files')
SILICO_FOLDER = os.path.join(HERE, '..', 'silico', 'test_silico_dicts')


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
def polyester_params():
    '''
    Parameters object for polyesters.

    Returns:
        Parameters: Parameters object.
    '''
    return generate_parameters(
        params_json=os.path.join(
            INPUTS_FOLDER,
            'full_input_parameters_polyester_LA_GA.json'))

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

@pytest.fixture
def test_polyester_silico_dict():
    '''
    Dict of l- and g-containing polyester sequences and their MS2 fragments,
    including signature ions.

    Returns:
        Dict[str, dict]: {sequence (str): MS2 fragments (dict)}
    '''
    polyester_file = os.path.join(
        SILICO_FOLDER,
        'polyester_l_g_silico_dict_ms2.json')
    with open(polyester_file) as f:
        load = json.load(f)
        return load

@pytest.mark.unit
def test_polyester_ms2(
    polyester_params,
    polyester_polymer,
    test_polyester_silico_dict
):
    '''
    Test for standard polyester c, i and j fragments.

    Args:
        polyester_params (Parameters): Parameters object.
        polyester_polymer (Polymer): Polymer object.
        test_polyester_silico_dict (Dict[str, dict]): MS2 fragment dict for
            l- and g-containing polyester sequences.
    '''
    #  make local copy of polyester params for changing parameters
    local_params = copy.deepcopy(polyester_params)

    #  ensure only sequences isomeric to one or more targets are generated
    local_params.silico.isomeric_targets = list(
        test_polyester_silico_dict.keys())

    #  generate ms2 dict from parameters
    ms2_dict = build_linear_fragments_sequence_dict(
        sequences=[seq for seq in test_polyester_silico_dict],
        params=polyester_params,
        polymer=polyester_polymer
    )

    #  iterate through test dict, and ensure that expected fragment dicts
    #  match between ms2_dict and test dict
    for seq, frags in test_polyester_silico_dict.items():
        for frag, info in frags.items():
            if frag != 'signatures':
                assert info == ms2_dict[seq][frag]

@pytest.mark.unit
def test_single_fragment_series(params, polymer, test_peptide_silico_dict):
    '''
    Tests linear fragment series for standard peptide sequences.

    Args:
        params (Parameters): parameters object.
        polymer (Polymer): polymer object.
        test_peptide_silico_dict (Dict[str, dict]): dict of standard peptide
            sequences and linear fragments.
    '''
    ms2_dict = build_linear_fragments_sequence_dict(
        sequences=[seq for seq in test_peptide_silico_dict],
        params=params,
        polymer=polymer
    )

    for seq, fragments in ms2_dict.items():
        for fragment, masses in fragments.items():
            assert masses == test_peptide_silico_dict[seq][fragment]


@pytest.mark.unit
def test_exception_fragments(
    depsi_params,
    depsi_polymer,
    test_depsi_silico_dict
):
    '''
    Tests linear fragment series for standard peptide sequences.

    Args:
        depsi_params (Parameters): parameters object.
        depsi_polymer (Polymer): polymer object.
        test_depsi_silico_dict (Dict[str, dict]): dict of depsipeptide
            sequences and linear fragments.
    '''
    ms2_dict = build_linear_fragments_sequence_dict(
        sequences=[seq for seq in test_depsi_silico_dict],
        params=depsi_params,
        polymer=depsi_polymer
    )

    for seq, fragments in ms2_dict.items():
        for fragment, masses in fragments.items():
            assert masses == test_depsi_silico_dict[seq][fragment]
