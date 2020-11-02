import pytest
from oligoss.postprocessing.postprocess_helpers import assign_confidence_sequence,\
    calculate_subsequence_coverage, assign_isomeric_groups, get_core_fragments

@pytest.fixture
def eg_confirmed_frags():
    return {"core": {
        "y3": [304.0792],
        "y4": [451.1476],
        "a1": [120.0813]},
        "signatures": {
            "ImF": [120.0813],
            "terminal_modifications": {
                "Ole": [69.1269]}}}

@pytest.fixture
def eg_silico_frags():
    return {
        "b1": [148.0763],
        "b2": [295.1447],
        "b3": [353.1502],
        "y2": [161.0087],
        "y3": [304.0792, 282.0978],
        "y4": [451.1476, 429.1662],
        "a1": [120.0813, 121.0892, 143.0705],
        "a2": [290.1389, 267.1497, 268.1576],
        "a3": [348.1444, 325.1552, 326.1631],
        "signatures": {
            "ImF": [120.0813],
            "terminal_modifications": {}}}

@pytest.mark.unit
def test_get_core_fragments():

    target_fragments = ['b1', 'b2', 'a1', 'a5', 'y1', 'y4', 'x3', 'x7', 'ImF']

    b_test = get_core_fragments(
        target_fragments=target_fragments, core_series=['b'])
    y_test = get_core_fragments(
        target_fragments=target_fragments, core_series=['y'])
    b_y_test = get_core_fragments(
        target_fragments=target_fragments, core_series=['b', 'y'])
    neg_test = get_core_fragments(
        target_fragments=target_fragments, core_series=['w'])

    assert sorted(b_test) == ['b1', 'b2']
    assert sorted(y_test) == ['y1', 'y4']
    assert sorted(b_y_test) == ['b1', 'b2', 'y1', 'y4']
    assert sorted(neg_test) == []

@pytest.mark.unit
def test_calculate_subsequence_coverage():

    single_block_test = calculate_subsequence_coverage(
        fragment_indices=[1, 2, 3, 4, 5, 7, 9, 11, 13])

    double_block_test = calculate_subsequence_coverage(
        fragment_indices=[1, 2, 3, 5, 7, 9, 10, 11, 12, 15, 19])

    negative_test = calculate_subsequence_coverage(
        fragment_indices=[1, 3, 5, 7, 9, 11, 13])

    assert single_block_test == 5
    assert double_block_test == 4
    assert negative_test == 1

@pytest.mark.unit
def test_assign_isomeric_groups():

    neg_test = assign_isomeric_groups(
        sequences=['ABCDE', 'ABCCC', 'ABDE', 'DDEE', 'BBAAB'])
    one_group_test = assign_isomeric_groups(
        sequences=['ABCDE', 'EDCBA', 'ACBDE', 'BDECA', 'ACDEB'])
    two_group_test = assign_isomeric_groups(
        sequences=['ABCDE', 'EDCBA', 'ACBBE', 'BBECA', 'ACBEB'])
    three_group_test = assign_isomeric_groups(
        sequences=['ABCDE', 'EDCBA', 'AABDE', 'BDEAA', 'ABDEB'])

    assert max(neg_test.values()) == 5
    assert max(one_group_test.values()) == 1
    assert max(two_group_test.values()) == 2
    assert max(three_group_test.values()) == 3

@pytest.mark.unit
def test_assign_confidence_sequence(eg_confirmed_frags, eg_silico_frags):

    b_y_positive_test_1 = assign_confidence_sequence(
        insilico_fragments=eg_silico_frags,
        confirmed_fragments=eg_confirmed_frags,
        core_fragment_series=['b', 'y'],
        optional_fragments=['b1'],
        exclude_fragments=None,
        essential_fragments=None,
        essential_fragment_cap=0,
        sequence_coverage_weight=0)

    b_y_positive_test_2 = assign_confidence_sequence(
        insilico_fragments=eg_silico_frags,
        confirmed_fragments=eg_confirmed_frags,
        core_fragment_series=['b', 'y'],
        optional_fragments=['b1'],
        exclude_fragments=None,
        essential_fragments=None,
        essential_fragment_cap=0,
        sequence_coverage_weight=1)

    b_y_a_positive_test = assign_confidence_sequence(
        insilico_fragments=eg_silico_frags,
        confirmed_fragments=eg_confirmed_frags,
        core_fragment_series=['b', 'y', 'a'],
        optional_fragments=['a1'],
        exclude_fragments=None,
        essential_fragments=None,
        essential_fragment_cap=0,
        sequence_coverage_weight=0)

    negative_test = assign_confidence_sequence(
        insilico_fragments=eg_silico_frags,
        confirmed_fragments=eg_confirmed_frags,
        core_fragment_series=['b'],
        optional_fragments=None,
        exclude_fragments=['b3'],
        essential_fragments=None,
        essential_fragment_cap=0,
        sequence_coverage_weight=0)

    assert b_y_positive_test_1 == 40.00
    assert b_y_positive_test_2 == 33.3333
    assert b_y_a_positive_test == 33.33
    assert negative_test == 0
