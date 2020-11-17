"""
Script to create specifically formatted files from the original probes file
Picard and merged,annotated BED file
"""

import subprocess
import sys
import logging
import os
import pandas as pd
from natsort import natsorted
from ngs_capture_qc.utils import check_probe_format
if sys.version_info[0] < 3: 
    from StringIO import StringIO
else:
    from io import StringIO

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('probefile', help='The probe file from the vendor file')
    parser.add_argument('refgene_bed', help="UCSC RefGene gene data in bed format, chrm|start|stop|gene")
    parser.add_argument('bedtools', default='',help='Path to bedtools, accepts binary or singularity image')
    parser.add_argument('outdir', default='.', help="Output directory for summary scripts")

def write_merged_bed(probes, bedtools, temp_merged_bed):
    """Given correctly formatted probes, write the merged bed file"""
    #First, create merged bed file
    write_probes=open(temp_merged_bed, 'w')
    merge_probes_args = [x for x in bedtools.split(' ')]+['bedtools', 'merge', '-i', probes]
    subprocess.call(merge_probes_args, stdout=write_probes) 
    write_probes.close()

def write_annotated_bed(temp_merged_bed, bedtools, refgene, anno_bed):
    """Given merged bed file, replacing the annotation with gene names"""
    print("running:",  temp_merged_bed, bedtools, refgene, anno_bed)
    #Next, annotate this file
    intersect_args = [x for x in bedtools.split(' ')]+['bedtools', 'intersect', '-a', temp_merged_bed, '-b', refgene, '-loj'] 
    annotate = subprocess.Popen(intersect_args, stdout=subprocess.PIPE)
    data = StringIO(annotate.communicate()[0].decode('utf-8'))
    df = pd.read_csv(data, sep='\t', header=None)
    #Ignore all columns in merged output except the 7 we care about
    df=df.iloc[:,:7]
    df.columns=['chrom','start','stop','Rchr','Rstart','Rstop','gene']
    #Drop duplicates based on chrm/start/stop of probes/bed, not of ucsc annotation (which is the Rchr/Rstart/Rstop)
    df.drop_duplicates(subset=['chrom','start','stop','gene'],inplace=True)
    #Now we need to drop duplicates that have different gene names but duplicate positions
#    df=df.groupby(['chrom','start','stop'])['gene'].apply('-'.join).reset_index()
    df=df.groupby(['chrom','start','stop']).gene.unique().apply(lambda x: ';'.join(x)).reset_index()

    df.replace(to_replace=r'^\.$', value='intergenic', regex=True, inplace=True)
    #Sort by chromosome and start
    #See https://stackoverflow.com/questions/29580978/naturally-sorting-pandas-dataframe
    df.chrom=df.chrom.astype('category')
    df.chrom.cat.reorder_categories(natsorted(set(df.chrom)), inplace=True, ordered=True)
    df.start = df.start.astype('category')
    df.start.cat.reorder_categories(natsorted(set(df.start)), inplace=True, ordered=True)
    df.sort_values('chrom', inplace=True)
    df.to_csv(anno_bed, columns=['chrom','start','stop','gene'],header=None,index=False, sep='\t')

def wrtie_antitarget_bed():
    """Given bed file, create bed file of regions not expected to be covered"""
    pass

def create_bed(probes,output_basename,refgene_bed, bedtools):
    """Inital step for new assay, validate probe file and write clean, annotated bed file"""
    #Write temp clean probe file for bedtools usage
    probes_temp=output_basename+'-TEMP.probes'
    probes.to_csv(probes_temp, columns=['chrom','start','stop'],header=False,sep='\t', index=False)

    #Write temp merged bed file
    temp_merged_bed=output_basename+'-TEMP.bed'
    write_merged_bed(probes_temp, bedtools, temp_merged_bed)

    #Write merged, annotated bed file
    anno_bed=output_basename+'.anno.bed'
    write_annotated_bed(temp_merged_bed, bedtools, refgene_bed, anno_bed)

    # #Remove temp merged bed file
    # os.remove(probes_temp)
    # os.remove(temp_merged_bed)

def create_picard_bed(probes, output_basename):
    """Use doc/PicardHeader and probe file to create a file 
    in the format required by picard
    """
    picard_bed=output_basename+'.Picard.bed'
    with open('docs/PicardHeader','r') as header:
        with open(picard_bed, 'w') as picard_out:
            for line in header:
                picard_out.write(line)
    probes.to_csv(picard_bed, columns=['chrom','start','stop','strand','annotation'],header=False,index=False,sep='\t', mode='a')

def action(args):
    #Read in the probes
    probes = pd.read_csv(open(args.probefile,'r'), delimiter='\t', header=None)

    #Assert the probes are in the correct format for processing
    probes = check_probe_format(probes)
    
    #setup name for resulting files
    probe_basename=os.path.splitext(os.path.basename(args.probefile))[0]
    output_basename=os.path.join(os.path.join(args.outdir,probe_basename))

    #Now, create files based on CLI arguments
    #Parse probes, write clean bed file
    if args.bedtools.endswith('img'):
        bedtools='singularity exec --bind {} --pwd {} {}'.format(os.getcwd(), os.getcwd(), args.bedtools)
    else:
        bedtools=args.bedtools
    create_bed(probes, output_basename, args.refgene_bed, bedtools)

    #Parse probes, write picard file
    create_picard_bed(probes, output_basename)
