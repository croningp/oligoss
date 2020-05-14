"""
This file contains tests for silico module helper functions. These are basic
functions required to perform simple tasks in the silico module, such as
calculate sequence masses or manipulate sequence strings.
"""

import os
import pytest

from ...utils.parameter_handlers.parameter_handlers import generate_parameters
from ...silico.polymer_info.polymer import Polymer
from ...silico.ms1_silico import generate_ms1_ions
from ...silico.helpers.helpers import (
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

    #  generate dict of sequences and corresponding MS1 precursor ions
    seq_ms1_dict = generate_ms1_ions(
        params=params,
        polymer=polymer,
        sequencing=True
    )

    #  generate dict of compositions and corresponding MS1 precursor ions
    composition_ms1_dict = generate_ms1_ions(
        params=params,
        polymer=polymer,
        sequencing=False
    )

    #  generate list of sequences
    sequences = generate_all_sequences(
        polymer=polymer,
        params=params,
        sequencing=True
    )

    #  check that sequence dict has all target sequences and that composition
    #  dict != sequence dict
    assert sorted(sequences) == sorted(seq_ms1_dict.keys())
    assert len(sequences) > len(composition_ms1_dict)

    #  iterate through sequences, making sure precursor m/z values are correct
    for sequence in sequences:
        precursors = seq_ms1_dict[sequence]
        composition = list(filter(
            lambda comp: sorted(comp) == sorted(sequence),
            composition_ms1_dict.keys()
        ))[0]
        comp_precursors = composition_ms1_dict[composition]

        calculated_precursors = ionize_sequence_precursors(
            sequence=sequence,
            params=params,
            polymer=polymer
        )

        assert precursors == comp_precursors == calculated_precursors
