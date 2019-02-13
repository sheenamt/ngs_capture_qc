"""
Test utils module.
"""

import logging

from ngs_capture_qc.utils import Opener

from __init__ import TestBase, get_testfile
log = logging.getLogger(__name__)


class TestOpener(TestBase):
    def setUp(self):
        with open(get_testfile('lorem.txt')) as f:
            self.firstline = next(f)

    def test01(self):
        for suffix in ['txt', 'gz', 'bz2']:
            fn = get_testfile('lorem.'+suffix)
            with Opener()(fn) as fobj:
                self.assertEqual(next(fobj), self.firstline)

    def test02(self):
        for suffix in ['txt', 'gz', 'bz2']:
            fn = get_testfile('lorem.'+suffix)
            with Opener('r')(fn) as fobj:
                self.assertEqual(next(fobj), self.firstline)
