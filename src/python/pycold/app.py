""" Main bootstrap and configuration module for pycold.  Any module that
requires configuration or services should import app and obtain the
configuration or service from here.
app.py enables a very basic but sufficient form of loose coupling
by setting names of services & configuration once and allowing other modules
that require these services/information to obtain them by name rather than
directly importing or instantiating.
Module level constructs are only evaluated once in a Python application's
lifecycle, usually at the time of first import. This pattern is borrowed
from Flask.
source: https://github.com/repository-preservation/lcmap-pyclass/blob/develop/pyclass/app.py
"""

import os
import yaml


class Defaults(dict):
    def __init__(self, config_path='parameters.yaml'):
        with open(config_path, 'r') as f:
            super(Defaults, self).__init__(yaml.safe_load(f.read()))

    def __getattr__(self, name):
        if name in self:
            return self[name]
        else:
            raise AttributeError('No such attribute: ' + name)

    def __setattr__(self, name, value):
        self[name] = value

    def __delattr__(self, name):
        if name in self:
            del self[name]
        else:
            raise AttributeError('No such attribute: ' + name)


defaults = Defaults(os.path.join(os.path.dirname(__file__), 'constants.yaml'))
