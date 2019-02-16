"""
Package description (in ngs_capture_qc/__init__.py)

Update me!
"""

import glob
from os import path

_data = path.join(path.dirname(__file__), 'data')


def package_data(fname, pattern=None):
    """Return the absolute path to a file included in package data,
    raising ValueError if no such file exists. If `pattern` is
    provided, return a list of matching files in package data
    (ignoring `fname`).

    """

    if pattern:
        return glob.glob(path.join(_data, pattern))

    pth = path.join(_data, fname)

    if not path.exists(pth):
        raise ValueError('Package data does not contain the file %s' % fname)

    return pth


try:
    with open(path.join(path.dirname(__file__), 'data', 'ver')) as f:
        __version__ = f.read().strip().replace('-', '+', 1).replace('-', '.')
except Exception as e:
    __version__ = ''
