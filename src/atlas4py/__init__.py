# -*- coding: utf-8 -*-

"""
Python bindings for the ECMWF Atlas mesh library.
"""
import atexit

from ._atlas4py import *
from ._atlas4py import __version__

_atlas4py._initialise()
atexit.register(lambda: _atlas4py._finalise())
