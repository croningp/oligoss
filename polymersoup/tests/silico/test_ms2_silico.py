"""
This file contains tests for silico module helper functions. These are basic
functions required to perform simple tasks in the silico module, such as
calculate sequence masses or manipulate sequence strings.
"""

import os
import pytest
import random

from ...utils.parameter_handlers.parameter_handlers import generate_parameters
from ...silico.polymer_info.polymer import Polymer
from ...silico.ms2_silico import (
    build_fragment_dict_for_sequences,
    build_linear_fragments_sequence_dict
)
from ...silico.helpers.helpers import generate_all_sequences

HERE = os.path.abspath(os.path.dirname(__file__))
INPUTS_FOLDER = os.path.join(HERE, '..', 'input_param_files')

@pytest.fixture
def params():
    return generate_parameters(
        params_json=os.path.join(
            INPUTS_FOLDER,
            'full_input_parameters_large_sequences.json'))

@pytest.fixture
def polymer(params):
    return Polymer(params_obj=params)

@pytest.fixture
def test_sequences(params, polymer):
    return generate_all_sequences(
        polymer=polymer,
        params=params,
        sequencing=True
    )

@pytest.mark.unit
def test_single_fragment_series(params, polymer, test_sequences):

    for fragment_series in polymer.fragment_info:

        fragment_dict = build_fragment_dict_for_sequences(
            sequences=random.sample(test_sequences, k=1000),
            params=params,
            polymer=polymer,
            fragment_series=fragment_series
        )

        assert fragment_dict

@pytest.mark.unit
def test_linear_fragment_dict_builder(params, polymer, test_sequences):

    linear_fragment_dict = build_linear_fragments_sequence_dict(
        sequences=random.sample(test_sequences, k=500),
        params=params,
        polymer=polymer
    )

    assert linear_fragment_dict
