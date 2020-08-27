import pytest

@pytest.fixture
def example_MS1_dict():
    """
    Example MS1 ripper data retrieved from MongoDB.

    Returns:
        List[Dict[str, Any]]: ripper MS1 dict.
    """
    return [
        {"spectrum_id": "spectrum_1", "spectrum": {
            "159.9641": 1843, "159.9644": 2790, "159.9646": 3108,
            "159.9649": 2129, "retention_time": 1.234}},
        {"spectrum_id": "spectrum_2", "spectrum": {
            "217.9568": 3373, "217.9573": 1906, "217.9586": 2677,
            "217.9590": 4325, "217.9594": 5098, "217.9599": 3508,
            "retention_time": 1.888}},
        {"spectrum_id": "spectrum_3", "spectrum": {
            "512.4443": 4955, "513.3324": 43495, "514.2045": 5023,
            "516.2305": 304395,
            "retention_time": 10.111}}
    ]

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
