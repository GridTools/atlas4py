# -*- coding: utf-8 -*-

"""
Python bindings for the ECMWF Atlas mesh library.
"""
import atexit

from ._version import __version__
from ._atlas4py import *

_atlas4py._initialise()
atexit.register(lambda: _atlas4py._finalise())
