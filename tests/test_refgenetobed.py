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

    def testGetIntList(self):
        parsed_list=refgene_to_bed.get_int_list(self.int_string)
        self.assertListEqual(parsed_list,self.int_list)

    def testGetStringList(self):
        parsed_string=refgene_to_bed.get_string_list(self.int_list)
        self.assertEqual(parsed_string,self.int_string)

 
    def testBedKey(self):
        indata={'chrom': '1', 'chromStart': '66999638', 'chromEnd': '67216822', 'name': 'SGIP1', 'refgene': 'NM_032291', 'exonCount': '25', 'exonSizes': '413,64,25,72,57,55,176,12,12,25,52,86,93,75,501,128,127,60,112,156,133,203,65,165,8067', 'exonStarts': '66999638,67091529,67098752,67101626,67105459,67108492,67109226,67126195,67133212,67136677,67137626,67138963,67142686,67145360,67147551,67154830,67155872,67161116,67184976,67194946,67199430,67205017,67206340,67206954,67208755', 'exonEnds': '67000051,67091593,67098777,67101698,67105516,67108547,67109402,67126207,67133224,67136702,67137678,67139049,67142779,67145435,67148052,67154958,67155999,67161176,67185088,67195102,67199563,67205220,67206405,67207119,67216822'}
        expected=['1', 66999638, 67216822]
        self.assertListEqual(refgene_to_bed.bed_key(indata), expected)

    def testRefGeneToBed(self):
        refgene=os.path.join(config.datadir, 'test.refGene')
        expected_output=os.path.join(config.datadir, 'expected-test.refGene.bed')
        test_output=os.path.join(self.outdir, 'test-output.refgene.bed')
        cmd=["/mnt/disk10/users/sheenams/ngs_capture_qc/cap_qc.py", "refgene_to_bed", refgene, test_output]
        subprocess.call(cmd)
        self.assertTrue(filecmp.cmp(expected_output, test_output))
