"""
Script to create specifically formatted files from the original probes file
"""

import subprocess
import sys
import logging
import shutil
import os
import pandas as pd
import numpy as np
from ngs_capture_qc.utils import chromosomes, check_probe_format
if sys.version_info[0] < 3: 
    from StringIO import StringIO
else:
    from io import StringIO

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('probefile', help='The probe file from the vendor file')
    parser.add_argument('-p','--picard',action='store_true',help="Create picard formatted file")
    parser.add_argument('-b','--bed', action='store_true', help="Create bed formatted file, bedtools also required")
    parser.add_argument('-g','--refseq_bed', help="UCSC RefSeq gene data in bed format, chrm|start|stop|gene")
    parser.add_argument('--bedtools', default='',help='Path to bedtools, accepts binary or singularity image')

def write_merged_bed(probes, bedtools, temp_merged_bed):
    """Given correctly formatted probes, write the merged bed file"""
    #First, create merged bed file
    write_probes=open(temp_merged_bed, 'w')
    merge_probes_args = [x for x in bedtools.split(' ')]+['bedtools', 'merge', '-i', probes]
    merge_probes = subprocess.call(merge_probes_args, stdout=write_probes) 
    write_probes.close()

def write_annotated_bed(temp_merged_bed, bedtools, refseq, anno_bed):
    """Given merged bed file, replacing the annotation with gene names"""
    #Next, annotate this file
    intersect_args = [x for x in bedtools.split(' ')]+['bedtools', 'intersect', '-a', temp_merged_bed, '-b', refseq, '-wa', '-wb'] 
    annotate = subprocess.Popen(intersect_args, stdout=subprocess.PIPE)
    data = StringIO(annotate.communicate()[0].decode('utf-8'))
    df = pd.read_csv(data, sep='\t', header=None)
    #Ignore all columns in merged output except the 7 we care about
    df=df.iloc[:,:7]
    df.columns=['chrom','start','stop','Rchr','Rstart','Rstop','gene']
    df.drop_duplicates(subset=['chrom','start','stop','gene'],inplace=True)
    df.to_csv(anno_bed, columns=['chrom','start','stop','gene'],header=None,index=False, sep='\t')

def create_bed(probes,probe_basename,refseq_bed, bedtools):
    """Inital step for new assay, validate probe file and write clean, annotated bed file"""
    #Write temp clean probe file for bedtools usage
    probes_temp=os.path.basename(probe_basename)+'-TEMP.probes'
    probes.to_csv(probes_temp, columns=['chrom','start','stop'],header=False,sep='\t', index=False)

    #Write temp merged bed file
    temp_merged_bed=os.path.basename(probe_basename)+'-TEMP.bed'
    write_merged_bed(probes_temp, bedtools, temp_merged_bed)

    #Write merged, annotated bed file
    anno_bed=os.path.basename(probe_basename)+'-anno.bed'
    write_annotated_bed(temp_merged_bed, bedtools, refseq_bed, anno_bed)
    
    #Remove temp merged bed file
    os.remove(probes_temp)
    os.remove(temp_merged_bed)

def create_picard_bed(probes, probe_basename):
    """Use doc/PicardHeader and probe file to create a file 
    in the format required by picard
    """
    picard_bed=os.path.basename(probe_basename)+'.Picard.bed'
    with open('docs/PicardHeader','r') as header:
        with open(picard_bed, 'w') as picard_out:
            for line in header:
                picard_out.write(line)
    probes.to_csv(picard_bed, columns=['chrom','start','stop','strand','annotation'],header=False,index=False,sep='\t', mode='a')

def action(args):
    #Read in the probes
    probes = pd.read_csv(open(args.probefile,'r'), delimiter='\t')

    #Assert the probes are in the correct format for processing
    probes = check_probe_format(probes)

    #setup name for resulting files
    probe_basename=os.path.splitext(os.path.basename(args.probefile))[0]

    #Now, create files based on CLI arguments
    #Parse probes, write clean bed file
    if args.bed:
        if args.refseq_bed and args.bedtools:
            if args.bedtools.endswith('img'):
                bedtools='singularity exec --bind {} --pwd {} {}'.format(os.getcwd(), os.getcwd(), args.bedtools)
            else:
                bedtools=args.bedtools
            create_bed(probes, probe_basename, args.refseq_bed, bedtools)
        else:
            raise Exception('ERROR:Creating the bed file requires a refseq gene file and the path/image of bedtools')

    #Parse probes, write picard file
    if args.picard:
        create_picard_bed(probes, probe_basename)

    #Create the filtered refseq file, based on prefered transcripts
