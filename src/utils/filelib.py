'''For checking the status of and creating files'''

from typing import Optional, Container
from pathlib import Path

from polymerist.genutils.decorators.functional import allow_string_paths


def validate_file_path(
        path : Path,
        check_missing : bool=False,
        check_already_exists : bool=False,
        check_has_extension : bool=False,
        valid_extensions : Optional[Container[str]]=None,
    ) -> None:
    '''Check that a path exists, is pathlike, and has a valid extension
    Performs no check by default if no arguments other than the path are passed; these NEED to be supplied when called!'''
    # Meta-error for nonsensical argument checks
    if check_missing and check_already_exists:
        raise ValueError('Incongruent checks requested; invalid to check both that a file already exists AND is missing')
    
    # Conditional file checks
    if check_missing and not path.exists():
        raise FileNotFoundError(f'No file exists at "{path}"')
    if check_already_exists and path.exists():
        raise PermissionError(f'File already exists at "{path}"')
    if check_has_extension and not path.suffix:
        raise IsADirectoryError(f'Input file missing file extension, appears like directory ("{path}")')
    if valid_extensions and (path.suffix not in valid_extensions):
        raise ValueError(
            f'Cannot read data from {path.suffix} file, choose one of' \
            f'the following valid extensions: {[ext for ext in valid_extensions]}' # NOTE: using comprehension instead of list() to give expected output for dicts
        )
