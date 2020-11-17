"""
Test the filter_refgene script
"""

import subprocess
import filecmp
import logging
import os
import pandas as pd
from ngs_capture_qc.subcommands import filter_refgene

from __init__ import TestBase
import __init__ as config
log = logging.getLogger(__name__)


class TestFilterRefGene(TestBase):
    """
    Test the create files script, which writes files based on probes and 
    UCSC refGene table info. 
    """
    def setUp(self):
        self.outdir = self.mkoutdir()

    def testReadRefGene(self):
        """Test the mapping of field headers to file"""
        refgene=open(os.path.join(config.datadir,'test.refGene'),'r')
        data=filter_refgene.read_refgene(refgene)
        test_data=[x for x in data if x['name']=='FOXA1']
        expected_data=dict([('bin', '875'), ('refgene', 'NM_004496'), ('chrom', '14'), ('strand', '-'), ('txStart', '38058756'), ('txEnd', '38064325'), ('cdsStart', '38060569'), ('cdsEnd', '38064177'), ('exonCount', '2'), ('exonStarts', '38058756,38064105,'), ('exonEnds', '38061916,38064325,'), ('score', '0'), ('name', 'FOXA1'), ('cdsStartStat', 'cmpl'), ('cdsEndStat', 'cmpl'), ('exonFrames', '0,0,')])
        refgene.close()
        self.assertDictEqual(dict(test_data[0]), expected_data)
        
    def testCheckOverlapping(self):
        feature1=[('ARID1A', 27022521, 27108601), ('ANGPTL1', 178818669, 178840215), ('ABL2', 179068461, 179112224), ('AKT3', 243663020, 244006886)]
        feature2=[('NAB2', 57482676, 57489259), ('STAT6', 57489186, 57505196)]

        #No overlapps in first list, returns empty list
        self.assertFalse(filter_refgene.check_overlapping(feature1))
        #Overlap in second list, returns overlap
        self.assertRaises(ValueError, filter_refgene.check_overlapping, feature2)

    def testFilterRefGene(self):
        refgene=os.path.join(config.datadir, 'test.refGene')
        genes=os.path.join(config.datadir, 'test.genes_for_filter')
        test_output=os.path.join(self.outdir, 'test-output-filtered.refgene')

        filtered_output=os.path.join(config.datadir, 'expected-filtered.refGene')
        cmd=["/mnt/disk10/users/sheenams/ngs_capture_qc/capqc", "filter_refgene", refgene, genes, test_output]        
        subprocess.call(cmd)
        self.assertTrue(filecmp.cmp(filtered_output, test_output))
