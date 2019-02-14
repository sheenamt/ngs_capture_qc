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
        test_data=[x for x in data if x['name2']=='MEGF6'][0]
        expected_data=dict([('bin', '76'), ('name', 'NM_001409'), ('chrom', '1'), ('strand', '-'), ('txStart', '3404505'), ('txEnd', '3528059'), ('cdsStart', '3407091'), ('cdsEnd', '3527832'), ('exonCount', '37'), ('exonStarts', '3404505,3407475,3409202,3410334,3410559,3410934,3411176,3412453,3413218,3413551,3413796,3414934,3415261,3415701,3416151,3416359,3417196,3417529,3417740,3418359,3421771,3421985,3422671,3424358,3425121,3425638,3426433,3427346,3428113,3428569,3431113,3431965,3440687,3496388,3511901,3519029,3527701,'), ('exonEnds', '3407153,3407523,3409331,3410463,3410688,3411063,3411305,3412582,3413347,3413683,3413925,3415063,3415390,3415830,3416280,3416488,3417328,3417658,3417872,3418485,3421906,3422120,3422800,3424487,3425253,3425809,3426556,3427466,3428251,3428692,3431236,3432091,3440810,3496493,3512011,3519164,3528059,'), ('score', '0'), ('name2', 'MEGF6'), ('cdsStartStat', 'cmpl'), ('cdsEndStat', 'cmpl'), ('exonFrames', '1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,0,')])
        refgene.close()
        self.assertDictEqual(test_data, expected_data)
        
    def testCheckOverlapping(self):
        feature1=[('ARID1A', 27022521, 27108601), ('ANGPTL1', 178818669, 178840215), ('ABL2', 179068461, 179112224), ('AKT3', 243663020, 244006886)]
        feature2=[('NAB2', 57482676, 57489259), ('STAT6', 57489186, 57505196)]

        #No overlapps in first list, returns empty list
        self.assertFalse(filter_refgene.check_overlapping(feature1))
        #Overlap in second list, returns overlap
        self.assertRaises(ValueError, filter_refgene.check_overlapping, feature2)

    def testFilterRefGene(self):
        refgene=os.path.join(config.datadir, 'test.refGene')
        genes=os.path.join(config.datadir, 'test.genes')
        test_output=os.path.join(self.outdir, 'test-output-filtered.refgene')

        filtered_output=os.path.join(config.datadir, 'expected-filtered.refGene')
        cmd=["/mnt/disk10/users/sheenams/ngs_capture_qc/cap_qc.py", "filter_refgene", refgene, genes, test_output]        
        subprocess.call(cmd)
        self.assertTrue(filecmp.cmp(filtered_output, test_output))
