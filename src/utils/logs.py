'''Configuration and contexts for job-to-file logging'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

import logging
LOGLEVELS_BY_NAME = logging.getLevelNamesMapping() # define once in header to avoid repreated calls

from typing import Any, Generator, Optional, Sequence, Union

import re
from pathlib import Path
from contextlib import contextmanager


# FORMATTING
TIMESTAMP_LOG = '%Y-%m-%d %H:%M:%S' # timestamp format to use for logging
LOG_FORMATTER = logging.Formatter('%(asctime)s.%(msecs)03d [%(levelname)-8s:%(module)16s:line %(lineno)-4d] - %(message)s', datefmt=TIMESTAMP_LOG) 

def format_error_for_log(error : Exception) -> str:
    '''Converts a raised Exception to a loggable string'''
    return f'{type(error).__name__}: {error!s}'

# FILE LOGGING
def create_file_logger(
    logfile_path : Union[str, Path],
    logger_name : str,
    level : int=logging.INFO,
    mode : str='a',
    formatter : logging.Formatter=LOG_FORMATTER,
    suppress_console_logs : bool=True,
) -> tuple[logging.Logger, logging.FileHandler]:
    '''Create a unique logger (bound to a file) for a Signac job'''

    # Create "mouthpiece" proxy logger which will handle log input
    logger = logging.getLogger(logger_name)
    if suppress_console_logs:
        ... # TODO: figure out how to suppress console logs (contextlib.redirect_stdout doesn't work!)

    # Create new handler to target file idempotently, returning a prior equivalent handler if one already exists
    for file_handler in logger.handlers:
        if file_handler.baseFilename == logfile_path:
            file_handler.mode = mode # TOSELF: not sure if this is safe?
            break # stop looking once a handler
    else:
        file_handler = logging.FileHandler(logfile_path, mode=mode) # make new handler if prior handler is not found
        
    # Configure handler properties
    file_handler.setLevel(level)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    return logger, file_handler

@contextmanager
def redirect_to_logfile(
        logfile_path : Union[str, Path],
        logger_name : str,
        level : int=logging.INFO,
        mode : str='a',
        formatter : logging.Formatter=LOG_FORMATTER,
        suppress_console_logs : bool=True,
        aux_loggers : Optional[Sequence[logging.Logger]]=None,
    ) -> Generator[logging.Logger, None, None]:
    '''Initializes a to-file logger for a job, and collects logs '''
    # Create a file handler
    proxy_logger, file_handler = create_file_logger(
        logfile_path=logfile_path,
        logger_name=logger_name,
        level=level,
        mode=mode,
        formatter=formatter,
        suppress_console_logs=suppress_console_logs,
    ) 

    # Temporarily bind handler to auxiliary loggers and suppress stdout
    if aux_loggers is None:
        aux_loggers = []

    orig_levels : dict[str, int] = {}
    for aux_logger in aux_loggers:
        orig_levels[aux_logger.name] = aux_logger.getEffectiveLevel() # record initial level for reversion at end
        if suppress_console_logs:
            ... # TODO: figure out how to suppress console logs (contextlib.redirect_stdout doesn't work!)
        aux_logger.addHandler(file_handler)

    # Execute 
    try:
        yield proxy_logger # supply mouthpiece for wrapped context to log to
    except Exception as e:
        proxy_logger.error(format_error_for_log(e))
    finally:
        for aux_logger in aux_loggers:
            aux_logger.setLevel(orig_levels[aux_logger.name])
            aux_logger.removeHandler(file_handler)

# LOG LEVEL SETTINGS
INT_REGEX = re.compile(r'^[-+]?\d+$') # define once in header to avoid repeated calls
def representable_as_int(inp : Any) -> bool: # NOTE: many, MANY ways to implement this; opted for RegEx, as it offers the most control over edge cases
    '''
    Check if an objects representation can be interpreted as an integer
    Integer as defined here consist of an optional +/- followed by at least one digit, and nothing else
    '''
    # try:
    #     _ = int(string) # overly inclusive to floats etc. which can be cast as ints
    #     return True
    # except ValueError:
    #     return False
    return (re.match(INT_REGEX, str(inp)) is not None) # more readable and (likely) quicker to skip isisntance pre-checks for str/int
    
def log_level_from_str(level_input : Union[int, str]) -> int:
    '''
    Parse a wide range of strings representing Python logging
    log levels into integers representing valid log levels
    '''
    if representable_as_int(level_input):
        level_input = int(level_input) # perform int conversion, if valid

    if isinstance(level_input, int):
        if level_input not in LOGLEVELS_BY_NAME.values():
            raise ValueError(f'Invalid logging level value {level_input}')
        return level_input
    elif isinstance(level_input, str):
        level_input = level_input.strip().upper() # permits level name inputs with leading/trailing whitespace and/or mixed letter case
        if (level_value := LOGLEVELS_BY_NAME.get(level_input, None)) is None:
            raise ValueError(f'No loggging level associated with name "{level_input}')
        return level_value
    else:
        raise TypeError(f'Cannot interpret object {level_input!r} of type {type(level_input).__name__} as log level')
loglevel_from_str = log_level_from_str # alias for convenience