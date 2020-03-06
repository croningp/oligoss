import pytest
import time
import json
from polymersoup.extractors.sequence_screening import (
    generate_EIC,
    generate_EIC_deprecated
)

@pytest.mark.unit
def test_generate_eic_optimization():
    print('Loading data...')
    RIPPER_DICT = '/home/group/data/polymermasssoup_data/ECK2_28_2_2_2.json'
    with open(RIPPER_DICT) as fd:
        ripper_dict = json.load(fd)
    args = {
        'ions': [458.1636, 475.1901, 480.145, 497.1715],
        'ms_level': 1,
        'err': 0.01,
        'err_abs': True,
        'ripper_dict': ripper_dict,
        'min_max_intensity': 0,
    }
    print('Running generate_EIC...\n')
    start_t = time.time()
    new_res = generate_EIC(**args)
    end_t = time.time()
    print(f'New time: {end_t - start_t:.2f} s')

    start_t = time.time()
    old_res = generate_EIC_deprecated(**args)
    end_t = time.time()
    print(f'Old time: {end_t - start_t:.2f} s')

    if not new_res:
        print('\nReturned list is empty []')

    assert old_res == new_res