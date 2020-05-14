import pytest
from ...extractors.data_extraction import match_mass

@pytest.fixture
def example_mass_list():
    return {"mass_list": [
            57.1412,
            124.4417,
            722.9747]}

@pytest.mark.unit
def test_mass_match(example_mass_list):
    exact_positive_test = match_mass(
        spectrum=example_mass_list,
        mass_range=[57.1412, 57.1412])

    positive_lower_bound_test = match_mass(
        spectrum=example_mass_list,
        mass_range=[57.1412, 57.1512])

    positive_upper_bound_test = match_mass(
        spectrum=example_mass_list,
        mass_range=[57.1312, 57.1412])

    negative_test = match_mass(
        spectrum=example_mass_list,
        mass_range=[125.4417, 125.5580])

    negative_lower_bound_test = match_mass(
        spectrum=example_mass_list,
        mass_range=[125.4316, 125.4416])

    negative_upper_bound_test = match_mass(
        spectrum=example_mass_list,
        mass_range=[124.4418, 124.4518])

    assert exact_positive_test == ['57.1412']
    assert positive_lower_bound_test == ['57.1412']
    assert positive_upper_bound_test == ['57.1412']
    assert negative_test == []
    assert negative_lower_bound_test == []
    assert negative_upper_bound_test == []
