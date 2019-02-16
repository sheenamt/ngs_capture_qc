"""
Test the summarize_assay script
"""

import subprocess
import filecmp
import logging
import os
import pandas as pd
from ngs_capture_qc.subcommands import summarize_assay

from __init__ import TestBase
import __init__ as config
log = logging.getLogger(__name__)


testfiles = os.path.join(config.datadir)


class TestSummarizeAssay(TestBase):
    """
    Test the create files script, which writes files based on probes and 
    UCSC refGene table info. 
    """

    def setUp(self):
        self.outdir = self.mkoutdir()
        self.assay=os.path.join(testfiles,'expected-ANNO.bed')

    def testCalculateTotalCovered(self):
        """Test calculation of probe coverage"""
        self.assertEqual(1504, summarize_assay.calculate_total_covered(self.assay))

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
        cmd=["/mnt/disk10/users/sheenams/ngs_capture_qc/cap_qc.py", "summarize_assay", self.assay, pref_trans, refgene, bedtools,"--outdir", self.outdir]
        subprocess.call(cmd)
        
        #files made in this include "overall_summary.txt", "per_refgene_summary.txt", "merged_probes.bed"
        expected_overall=os.path.join(testfiles, "expected-overall_summary.txt")
        expected_per_refgene=os.path.join(testfiles, "expected-per_refgene_summary.txt")
        self.assertTrue(filecmp.cmp(expected_overall, os.path.join(self.outdir, "overall_summary.txt")))
        self.assertTrue(filecmp.cmp(expected_per_refgene, os.path.join(self.outdir, "per_refgene_summary.txt")))




