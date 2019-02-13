"""
Test the summarize_assay script
"""

import subprocess
import filecmp
import logging
import os
import pandas as pd
from ngs_capture_qc.subcommands import create_files

from __init__ import TestBase
import __init__ as config
log = logging.getLogger(__name__)


assay_testfiles = os.path.join(config.datadir, 'assay_files')


class TestCreateFiles(TestBase):
    """
    Test the create files script, which writes files based on probes and 
    UCSC refGene table info. 
    """

    def setUp(self):
        self.outdir = self.mkoutdir()
        self.data=[['1', '3407032', '3407152','275744_14961323_MEGF6_chr1:3407091-3407153_1', '+'],['1', '3407036', '3407158','chr1:3407091-3407153_1', '-']]

    def break_data(self,data,string_to_break,index):
        self.new_data=[]
        for x in data:
            x[index]=string_to_break
            self.new_data.append(x)
        return self.new_data


    def testCheckFormat1(self):
        """Check that there is no header"""
        header_probes=pd.DataFrame(data=self.data,columns=['chrm','start','stop','annot','strand'])
        self.assertRaises(ValueError,create_files.check_format, header_probes)

    def testCheckFormat2(self):
        """Check that there are 5 columns """
        len_probes=pd.DataFrame(data=[x[:-1] for x in self.data])
        self.assertRaises(AssertionError,create_files.check_format,len_probes)

    def testCheckFormat3(self):
        """Check that the start and stop columns are integers """

        start_probes=pd.DataFrame(data=self.break_data(self.data,'bad',1))
        self.assertRaises(ValueError,create_files.check_format,start_probes)

        stop_probes=pd.DataFrame(data=self.break_data(self.data,'bad',2))
        self.assertRaises(ValueError, create_files.check_format,stop_probes)

    def testCheckFormat4(self):
        """Check that there is a strand column """

        strand_probes=pd.DataFrame(data=self.break_data(self.data,'bad',4))
        self.assertRaises(ValueError, create_files.check_format,strand_probes)

    def testCheckFormat5(self):
        """Check that there is a string of some sort in the annotation column """

        annot_probes=pd.DataFrame(data=self.break_data(self.data,'4',4))
        self.assertRaises(ValueError, create_files.check_format, annot_probes)

            
