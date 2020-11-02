"""
This file contains tests for ms1 silico module functions. These are functions
used to generate precursor ions. NOTE: all functions tested here do NOT
interface with MongoDB via pymongo. There is a separate _db file for ms1 silico
functions that use pymongo.
"""

import os
import pytest

from oligoss.utils.parameter_handlers.parameter_handlers import generate_parameters
from oligoss.silico.polymer_info.polymer import Polymer
from oligoss.silico.ms1_silico import generate_ms1_ions
from oligoss.silico.helpers.helpers import (
    generate_all_sequences,
    ionize_sequence_precursors
)

HERE = os.path.abspath(os.path.dirname(__file__))
INPUTS_FOLDER = os.path.join(HERE, '..', 'input_param_files')

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
def test_generate_ms1_ions(params, polymer):
    """
    Tests MS1 precursor dicts.

    Args:
        params (Parameters): Parameters object
        polymer (Polymer): Polymer object
    """

    #  create dict of sequences and corresponding MS1 precursor ions
    seq_ms1 = generate_ms1_ions(
        params=params,
        polymer=polymer,
        sequencing=True
    )
    seq_ms1 = {seq: precursor for (seq, precursor) in seq_ms1}

    #  generator of compositions and corresponding MS1 precursor ions
    composition_ms1 = generate_ms1_ions(
        params=params,
        polymer=polymer,
        sequencing=False
    )
    composition_ms1 = {seq: precursor for (seq, precursor) in composition_ms1}

    #  generate list of sequences
    sequences = generate_all_sequences(
        polymer=polymer,
        params=params,
        sequencing=True
    )
    sequences = sorted([seq for seq in sequences])

    #  check that sequence dict has all target sequences and that composition
    #  dict != sequence dict
    assert sequences == sorted(seq_ms1.keys())
    assert len(sequences) > len(composition_ms1)

    #  check that there are no repeat sequences
    assert sorted(list(set(sequences))) == sequences

    #  iterate through sequences, making sure precursor m/z values are correct
    for seq in sequences:

        precursors = ionize_sequence_precursors(
            sequence=seq,
            params=params,
            polymer=polymer
        )

        composition = [
            comp for comp in composition_ms1 if sorted(seq) == sorted(comp)][0]
        assert (
            precursors == composition_ms1[composition] == seq_ms1[seq])
