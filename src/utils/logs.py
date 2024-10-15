'''Configuration and contexts for job-to-file logging'''

import logging
from typing import Generator, Optional, Sequence, Union

from pathlib import Path
from contextlib import contextmanager


TIMESTAMP_LOG = '%Y-%m-%d %H:%M:%S' # timestamp format to use for logging
LOG_FORMATTER = logging.Formatter('%(asctime)s.%(msecs)03d [%(levelname)-8s:%(module)16s:line %(lineno)-4d] - %(message)s', datefmt=TIMESTAMP_LOG) 

def format_error_for_log(error : Exception) -> str:
    '''Converts a raised Exception to a loggable string'''
    return f'{type(error).__name__}: {error!s}'

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