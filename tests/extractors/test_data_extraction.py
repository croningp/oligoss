import pytest
from oligoss.extractors.data_extraction import retrieve_retention_times,\
    confirm_fragment

@pytest.fixture
def example_MS1_dict():
    return {"spectrum_1": {
        "159.9641": 1843, "159.9644": 2790, "159.9646": 3108, "159.9649": 2129,
        "retention_time": 1.234},
        "spectrum_2": {
            "217.9568": 3373, "217.9573": 1906, "217.9586": 2677,
            "217.9590": 4325, "217.9594": 5098, "217.9599": 3508,
            "retention_time": 1.888},
        "spectrum_3": {
            "512.4443": 4955, "513.3324": 43495, "514.2045": 5023,
            "516.2305": 304395,
            "retention_time": 10.111}}

@pytest.fixture
def example_MS2_dict():
    return {"spectrum_1": {
        "50.5803": 1915, "55.9342": 2547, "56.0492": 6106,
        "58.1675": 2080, "59.0132": 1777, "retention_time": "1.738227561067",
        "mass_list": [50.5803, 55.9342, 56.0492, 58.1675, 59.0132],
        "parent": "228.4263"},
        "spectrum_2": {
            "101.5577": 13530, "102.0591": 2931, "107.0600": 5300,
            "108.0554": 13671, "109.0630": 11450,
            "retention_time": "1.742648580783",
            "mass_list": [101.5577, 102.0591, 107.06, 108.0554, 109.063],
            "parent": "247.7614"},
        "spectrum_3": {
            "196.0708": 18939, "202.1081": 7813, "214.0814": 18958,
            "214.1077": 2617,
            "retention_time": "1.745956381317",
            "mass_list": [196.0708, 202.1081, 214.0814, 214.1077],
            "parent": "228.4263"}}

@pytest.mark.unit
def test_retrieve_retention_times(example_MS1_dict):
    test_rts = retrieve_retention_times(example_MS1_dict)
    assert test_rts == [1.234, 1.888, 10.111]

@pytest.mark.unit
def test_confirm_fragment(example_MS2_dict):
    pos_masses = (50.5803, 102.0591, 214.1077)

    abs_matches, positive_test_abs = confirm_fragment(
        masses=pos_masses,
        error=0.01,
        error_units='abs',
        spectra=example_MS2_dict)

    ppm_matches, positive_test_ppm = confirm_fragment(
        masses=pos_masses,
        error=1,
        error_units='ppm',
        spectra=example_MS2_dict)

    neg_matches, neg_test_abs = confirm_fragment(
        masses=[1000.00, 23.3040],
        error=0.01,
        error_units='abs',
        spectra=example_MS2_dict)

    assert list(positive_test_ppm.keys()) == [
        'spectrum_1', 'spectrum_2', 'spectrum_3']
    assert list(positive_test_abs.keys()) == [
        'spectrum_1', 'spectrum_2', 'spectrum_3']
    assert abs_matches, ppm_matches == pos_masses
    assert neg_matches == []
    assert neg_test_abs == {}
