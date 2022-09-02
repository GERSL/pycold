__version__ = '0.1.0'

# Ensure monkey patches run first
from . import monkey  # NOQA

from ._pycold_cython import *  # NOQA
