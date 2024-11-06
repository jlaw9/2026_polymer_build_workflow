'''Test robustness of check for integral representability'''

from typing import Any

import pytest
from ..logs import representable_as_int


repr_int_valid = [ # cases we expect to be representable as integers
    ## test integers, including null and signed values
    0,
    42,
    +42,
    -42,
    1000000,
    ## string analogues of the above ints
    '0',
    '42',
    '+42',
    '-42',
    '1000000',
] 

repr_int_invalid = [ ## cases we expect NOT to be representable as integers
    # float-like entries
    .0,
    0.,
    0.1,
    1.,
    .1,
    1.1,
    +1.1,
    -1.1,
    # string analogues of float-like entries
    '.0',
    '0.',
    '0.0',
    '1.',
    '.1',
    '1.1',
    '+1.1',
    '-1.1',
    # strings with no int analogue
    '',
    '-+1',
    '+-1',
    '12-34',
    '1.1.1',
    '1.1.0',
    '1.0.1',
    '1.0.0',
    '1.0.',
    '1..0',
    '1..',
    '1_000_000', # though valid in Python, will consider out-of-scope as input is tricky
    '1,000,000', # though valid typographically, will consider out-of-scope as input is tricky
    # textual strings
    'one',
    'forty-two',
    'bogus',
    '10words01',
    # overtly non-integer or string objects
    object(),
    (1,2,3),
    [1,2,3],
    {'one':'two'},
] 

@pytest.mark.parametrize('valid_input', repr_int_valid)
def test_int_reprs_valid(valid_input : Any) -> None:
    '''Test that valid int-representable objects are correctly recognized'''
    assert representable_as_int(valid_input) 

@pytest.mark.parametrize('invalid_input', repr_int_invalid)
def test_int_reprs_invalid(invalid_input : Any) -> None:
    '''Test that invalid int-representable objects are correctly rejected'''
    assert not representable_as_int(invalid_input)