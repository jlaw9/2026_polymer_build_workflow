'''Test logging level parsing'''

import logging
from typing import Union

import pytest

from ..logs import log_level_from_str


loglevel_expected_returns = {
    # integral values
    logging.CRITICAL : 50, # these are just shorthands for integers, but make clearer their role here
    logging.ERROR : 40,
    logging.WARN : 30,
    logging.INFO : 20,
    logging.DEBUG : 10,
    logging.NOTSET : 0,
    # integral values as strings
    '50' : 50,
    '40' : 40,
    '30' : 30,
    '20' : 20,
    '10' : 10,
    '0' : 0,
    # logging keywords
    'CRITICAL' : 50,
    'FATAL' : 50,
    'ERROR' : 40,
    'WARNING' : 30,
    'WARN' : 30,
    'INFO' : 20,
    'DEBUG' : 10,
    'NOTSET' : 0,
    # logging keywords lowercase
    'critical': 50,
    'fatal': 50,
    'error': 40,
    'warning': 30,
    'warn': 30,
    'info': 20,
    'debug': 10,
    'notset': 0,
    # logging keywords, mixed case + extra whitespace
    ' CRiTIcAl  ': 50,
    '  fATAl  ': 50,
    'eRror ': 40,
    '  warNINg': 30,
    '   WArN': 30,
    '  Info ': 20,
    'debUg ': 10,
    ' notSEt': 0,
}
loglevel_invalid = {
    # invalid integers/numerics
    55,
    0.0,
    314,
    -10,
    # invalid integers as strings
    '55',
    '0.0',
    '314',
    '-10',
    # deliberate typos
    'critcical',
    'ifno',
    'INF',
    'debugg ',
    'wanr',
    'wanring',
}

@pytest.mark.parametrize('level_input, expected_output', loglevel_expected_returns.items())
def test_loglevels_outputs(level_input : Union[str, int], expected_output : int) -> None:
    '''Check that valid inputs for logging levels produce the expected outputs'''
    assert log_level_from_str(level_input) == expected_output

@pytest.mark.xfail(raises=(TypeError, ValueError), strict=True)
@pytest.mark.parametrize('level_input', loglevel_invalid)
def test_loglevels_invalid(level_input : Union[str, int]) -> None:
    '''Check that invalid inputs for logging levels correctly raise Exception'''
    _ = log_level_from_str(level_input)
