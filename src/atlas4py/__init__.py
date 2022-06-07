# -*- coding: utf-8 -*-

"""
Python bindings for the ECMWF Atlas mesh library.
"""

from ._version import __version__
from ._atlas4py import *

from .pyvista import *

def make_view(field):
    import numpy as np
    return np.array(field, copy=False)

