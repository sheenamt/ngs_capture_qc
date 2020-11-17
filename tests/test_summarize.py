"""
Test the summarize_assay script
"""

import subprocess
import filecmp
import logging
import os
import pandas as pd
from ngs_capture_qc.subcommands import summarize_assay
from ngs_capture_qc.utils import mkdir
#from __init__ import TestCaseSuppressOutput, TestBase,
#from __init__ import datadir as datadir
import __init__ as config
log = logging.getLogger(__name__)

import unittest

testfiles = 'testfiles'
outputdir = 'test_output'
mkdir(outputdir)

class TestSummarizeAssay(unittest.TestCase): #TestBase):
    """
    Test the create files script, which writes files based on probes and 
    UCSC refGene table info. 
    """
    outputdir = outputdir
    def mkoutdir(self, clobber=True):
        """
        Create outdir as outpudir/module.class.method (destructively
        if clobber is True).
        """
    
        funcname = '.'.join(self.id().split('.')[-3:])
        outdir = os.path.join(self.outputdir, funcname)
        mkdir(outdir, clobber)
        return outdir

    def setUp(self):
        self.outdir = self.mkoutdir()
        self.assay=os.path.join(testfiles,'expected-ANNO.bed')

    def testCalculateTotalCovered(self):
        """Test calculation of probe coverage"""
        self.assertEqual(1504, summarize_assay.calculate_total_covered(self.assay))

    def testExonTracker1(self):
        """Test exon parsing when interval completely within exon"""
        ES=['100','300','500','700']
        EE=['200','400','600','800']
        exons=summarize_assay.exonTracker(ES, EE)
        exons.insert(150,175)
        expected_exons={('100', '200'): True, ('300', '400'): False, ('500', '600'): False, ('700', '800'): False}
        self.assertDictEqual(exons.exons, expected_exons)

    def testExonTracker2(self):
        """Test exon parsing when interval starts before exon and ends within"""
        ES=['100','300','500','700']
        EE=['200','400','600','800']
        exons=summarize_assay.exonTracker(ES, EE)
        exons.insert(250,350)
        expected_exons={('100', '200'): False, ('300', '400'): True, ('500', '600'): False, ('700', '800'): False}
        self.assertDictEqual(exons.exons, expected_exons)
        
    def testExonTracker3(self):
        """Test exon parsing when interval starts within exon and ends after"""
        ES=['100','300','500','700']
        EE=['200','400','600','800']
        exons=summarize_assay.exonTracker(ES, EE)
        exons.insert(550,650)
        expected_exons={('100', '200'): False, ('300', '400'): False, ('500', '600'): True, ('700', '800'): False}
        self.assertDictEqual(exons.exons, expected_exons)
                
    def testExonTracker4(self):
        """Test exon parsing when interval starts before exon and ends after"""
        ES=['100','300','500','700']
        EE=['200','400','600','800']
        exons=summarize_assay.exonTracker(ES, EE)
        exons.insert(675,850)
        expected_exons={('100', '200'): False, ('300', '400'): False, ('500', '600'): False, ('700', '800'): True}
        self.assertDictEqual(exons.exons, expected_exons)

    def testExonTracker5(self):
        """Test exon parsing when interval starts before exon and ends after"""
        ES=['100','300','500','700']
        EE=['200','400','600','800']
        exons=summarize_assay.exonTracker(ES, EE)
        exons.insert(150,175)
        exons.insert(250,350)
        exons.insert(675,850)
        exons.insert(550,650)
        expected_exons={('100', '200'): True, ('300', '400'): True, ('500', '600'): True, ('700', '800'): True}
        self.assertDictEqual(exons.exons, expected_exons)



    def testSummarizeAssay(self):
        #Test file includes:
        # Region that covers single exon, split in two intervals (RPL10)
        # Region that covers multiple exons, (FOXA1)
        # Region that covers multiple exons, split in two intervals (MEGF6)
        # incorrect refseq/gene name mapping, (FAKE)
        # Gene in preferred trans that is not in probes, (GPR146)
        # Intergenic probe (chr2:47617462-47617582)
        pref_trans = os.path.join(testfiles, 'test.genes_for_summarize')
        refgene = os.path.join(testfiles, 'test.refGene.bed')
        bedtools='/mnt/disk2/com/container-images/bedtools-2.26.img'
        cmd=["/mnt/disk10/users/sheenams/ngs_capture_qc/capqc", "summarize_assay", self.assay, pref_trans, refgene, bedtools,"--outdir", self.outdir]
        subprocess.call(cmd)
        
        #files made in this include "overall_summary.txt", "per_refgene_summary.txt", "merged_probes.bed"
        expected_overall=os.path.join(testfiles, "expected-overall_summary.txt")
        expected_per_refgene=os.path.join(testfiles, "expected-pref_refgene_summary.txt")
        expected_per_refgene=os.path.join(testfiles, "expected-other_refgene_summary.txt")
        self.assertTrue(filecmp.cmp(expected_overall, os.path.join(self.outdir, "overall_summary.txt")))
        self.assertTrue(filecmp.cmp(expected_pref_refgene, os.path.join(self.outdir, "preferred_refgene_summary.txt")))
        self.assertTrue(filecmp.cmp(expected_other_refgene, os.path.join(self.outdir, "other_refgene_summary.txt")))

