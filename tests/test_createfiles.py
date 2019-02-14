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

testfiles=config.datadir

class TestCreateFiles(TestBase):
    """
    Test the create files script, which writes files based on probes and 
    UCSC refGene table info. 
    """

    def setUp(self):
        self.outdir = self.mkoutdir()
        self.output_basename=os.path.join(testfiles,'testoutput')
        self.data=[['1', '3407032', '3407152','275744_14961323_MEGF6_chr1:3407091-3407153_1', '+'],['1', '3407036', '3407158','chr1:3407091-3407153_1', '-']]
        self.bedtools_image='/mnt/disk2/com/container-images/bedtools-2.26.img'
        self.bedtools='singularity exec --bind {} --pwd {} {}'.format(os.getcwd(), os.getcwd(), self.bedtools_image)
        self.refgene_bed = os.path.join(testfiles, 'expected.refGene.bed')
        self.probe_file=os.path.join(testfiles,'test.probes')
        self.probes_df = None
        self.probes_temp = None
        self.setupTempProbes()

    def setupTempProbes(self):
        pf=open(self.probe_file,'r')
        probes = pd.read_csv(pf, delimiter='\t')
        self.probes_df = create_files.check_probe_format(probes)
        self.probes_temp=os.path.join(self.outdir,'testoutput-TEMP.probes')
        self.probes_df.to_csv(self.probes_temp, columns=['chrom','start','stop'],header=False,sep='\t', index=False)
        pf.close()

    def break_data(self,data,string_to_break,index):
        self.new_data=[]
        for x in data:
            x[index]=string_to_break
            self.new_data.append(x)
        return self.new_data

    def testCheckFormat1(self):
        """Check that there is no header"""
        header_probes=pd.DataFrame(data=self.data,columns=['chrm','start','stop','annot','strand'])
        self.assertRaises(ValueError,create_files.check_probe_format, header_probes)

    def testCheckFormat2(self):
        """Check that there are 5 columns """
        len_probes=pd.DataFrame(data=[x[:-1] for x in self.data])
        self.assertRaises(AssertionError,create_files.check_probe_format,len_probes)

    def testCheckFormat3(self):
        """Check that the start and stop columns are integers """
        start_probes=pd.DataFrame(data=self.break_data(self.data,'bad',1))
        self.assertRaises(ValueError,create_files.check_probe_format,start_probes)

        stop_probes=pd.DataFrame(data=self.break_data(self.data,'bad',2))
        self.assertRaises(ValueError, create_files.check_probe_format,stop_probes)

    def testCheckFormat4(self):
        """Check that there is a strand column """
        strand_probes=pd.DataFrame(data=self.break_data(self.data,'bad',4))
        self.assertRaises(ValueError, create_files.check_probe_format,strand_probes)

    def testCheckFormat5(self):
        """Check that there is a string of some sort in the annotation column """
        annot_probes=pd.DataFrame(data=self.break_data(self.data,'4',4))
        self.assertRaises(ValueError, create_files.check_probe_format, annot_probes)

    def testCreatePicardBed(self):
        """Use doc/PicardHeader and probe file to create a file 
        in the format required by picard
        """
        expected_output=os.path.join(testfiles,'expected.Picard.bed')
        testing_output=os.path.join(self.outdir,'testoutput')
        create_files.create_picard_bed(self.probes_df, testing_output)
        self.assertTrue(filecmp.cmp(expected_output, os.path.join(self.outdir, 'testoutput.Picard.bed')))
            

    def testWriteMerged_Bed(self):
        expected_output=os.path.join(testfiles, 'expected-TEMP.bed')
        testing_output=os.path.join(self.outdir,'testoutput.bed')
        create_files.write_merged_bed(self.probes_temp, self.bedtools, testing_output)
        self.assertTrue(filecmp.cmp(expected_output, testing_output))

    def testWriteAnnotatedBed(self):
        #Write merged, annotated bed file
        expected_output=os.path.join(testfiles,'expected-ANNO.bed')
        input_bed=os.path.join(testfiles, 'expected-TEMP.bed')
        anno_bed=os.path.join(self.outdir,'testoutput.anno.bed')
        create_files.write_annotated_bed(input_bed, self.bedtools, self.refgene_bed, anno_bed)
        self.assertTrue(filecmp.cmp(expected_output, anno_bed))

    def testCreateFiles(self):
        #Test running of whole script
        cmd=["/mnt/disk10/users/sheenams/ngs_capture_qc/cap_qc.py", "create_files", self.probe_file, "--refgene_bed", self.refgene_bed, "--outdir", self.outdir, "--bedtools", self.bedtools_image, "--picard", "--bed"]
        subprocess.call(cmd)
        
        expected_picard=os.path.join(testfiles,'expected.Picard.bed')
        expected_bed=os.path.join(testfiles,'expected-ANNO.bed')
        testing_picard=os.path.join(self.outdir,'test.Picard.bed')
        testing_bed=os.path.join(self.outdir,'test.anno.bed')
        self.assertTrue(filecmp.cmp(expected_picard, testing_picard))
        self.assertTrue(filecmp.cmp(expected_bed, testing_bed))
  
  

