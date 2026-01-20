'''Source code for high-throughput polymer structure generation projects'''

__author__ = 'Timotej Bernat'
__email__ = 'timotej.bernat@colorado.edu'

from pathlib import Path
_parent_dir = Path(__file__).parent.resolve()

try:
    from .utils.forcefields import FF_DIRNAME
except ImportError:
    from utils.forcefields import FF_DIRNAME

FF_DIR = _parent_dir / FF_DIRNAME
FF_DIR.mkdir(exist_ok=True)