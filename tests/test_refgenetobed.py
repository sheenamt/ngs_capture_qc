"""
Test the refgene_to_bed script
"""

import subprocess
import filecmp
import logging
import os
import pandas as pd
from ngs_capture_qc.subcommands import refgene_to_bed

from __init__ import TestBase
import __init__ as config

log = logging.getLogger(__name__)

testfiles = config.datadir

class TestRefGeneToBed(TestBase):
    """
    Test the refgene to bed file, which parses the 
    UCSC refGene table info into bed format
    """
    def setUp(self):
        self.outdir = self.mkoutdir()
        self.int_string='66999638,67091529,67098752,67101626'
        self.int_list=[66999638, 67091529, 67098752, 67101626]

    def testRefGeneToBed(self):
        refgene=os.path.join(config.datadir, 'test.refGene')
        expected_output=os.path.join(config.datadir, 'expected.refGene.bed')
        test_output=os.path.join(self.outdir, 'test-output.refgene.bed')
        cmd=["/mnt/disk10/users/sheenams/ngs_capture_qc/cap_qc.py", "refgene_to_bed", refgene, test_output]
        subprocess.call(cmd)
        self.assertTrue(filecmp.cmp(expected_output, test_output))
