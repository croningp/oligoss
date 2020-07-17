"""
This file contains tests for silico module helper functions. These are basic
functions required to perform simple tasks in the silico module, such as
calculate sequence masses or manipulate sequence strings.
"""

import os
import copy
import json
import pytest

from ...utils.parameter_handlers.parameter_handlers import generate_parameters
from ...silico.polymer_info.polymer import Polymer
from ...silico.ms2_silico import build_linear_fragments_sequence_dict
from ...silico.silico_handler import get_ms2_silico_dict_from_compositions
from ...silico.helpers.helpers import generate_all_sequences

HERE = os.path.abspath(os.path.dirname(__file__))
INPUTS_FOLDER = os.path.join(HERE, '..', 'input_param_files')
SILICO_FOLDER = os.path.join(HERE, '..', 'silico', 'test_silico_dicts')

@pytest.fixture
def params():
    return generate_parameters(
        params_json=os.path.join(
            INPUTS_FOLDER,
            'full_input_parameters_large_sequences.json'))

@pytest.fixture
def depsi_params():
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
        sequencing=True
    )
@pytest.fixture
def test_full_ms2_depsi_dict():
    fp = os.path.join(SILICO_FOLDER, "full_peptide_MSMS_silico.json")
    with open(fp) as f:
        load = json.load(f)
        return load

@pytest.mark.unit
def test_full_fragment_dict_depsi(
    depsi_params,
    depsi_polymer,
    test_full_ms2_depsi_dict
):
    """
    Tests full MS2 fragment builder for linear depsipeptide sequences: a, b, y
    fragments + signature and composition tags.

    Args:
        depsi_params (Parameters): parameters object (fixture).
        depsi_polymer (Polymer): polymer object (fixture).
        test_full_ms2_depsi_dict (Dict[str, dict]): dict of depsipeptide
            sequences and full MS2 fragments (fixture).
    """
    local_params = copy.deepcopy(depsi_params)
    local_params.silico.isomeric_targets = [
        seq for seq in test_full_ms2_depsi_dict]

    ms2_dict = get_ms2_silico_dict_from_compositions(
        params=local_params,
        polymer=depsi_polymer,
        ms1_hits=[seq for seq in test_full_ms2_depsi_dict]
    )
    for seq, frags in ms2_dict.items():
        for frag, info in frags.items():
            assert info == test_full_ms2_depsi_dict[seq][frag]


@pytest.mark.unit
def test_single_fragment_series(params, polymer, test_peptide_silico_dict):
    """
    Tests linear fragment series for standard peptide sequences.

    Args:
        params (Parameters): parameters object.
        polymer (Polymer): polymer object.
        test_peptide_silico_dict (Dict[str, dict]): dict of standard peptide
            sequences and linear fragments.
    """
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
    """
    Tests linear fragment series for standard peptide sequences.

    Args:
        depsi_params (Parameters): parameters object.
        depsi_polymer (Polymer): polymer object.
        test_depsi_silico_dict (Dict[str, dict]): dict of depsipeptide
            sequences and linear fragments.
    """
    ms2_dict = build_linear_fragments_sequence_dict(
        sequences=[seq for seq in test_depsi_silico_dict],
        params=depsi_params,
        polymer=depsi_polymer
    )

    for seq, fragments in ms2_dict.items():
        for fragment, masses in fragments.items():
            assert masses == test_depsi_silico_dict[seq][fragment]
