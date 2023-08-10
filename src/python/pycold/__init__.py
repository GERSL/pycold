__version__ = '0.1.2'

# Ensure monkey patches run first
from . import monkey  # NOQA

# from . import _colds_cython  # NOQA

# from ._colds_cython import (  # NOQA
#     test_func, _cold_detect,
#     obcold_reconstruct, sccd_detect, sccd_update)

from . import colds
# from . import _colds_cython
from .colds import cold_detect, obcold_reconstruct, sccd_detect, sccd_update, sccd_identify, calculate_sccd_cm
