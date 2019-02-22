"""
Convert UCSC refgene.txt files to BED format, for use in summarize_assay script
"""
 
import os
import sys 
import csv
from operator import itemgetter

import pandas
from natsort import natsorted
from ngs_capture_qc.utils import chromosomes
 
def build_parser(parser):
    parser.add_argument('refgene', help='UCSC table browser download')
    parser.add_argument('outfile', help='Output file', default = sys.stdout)

 
def action(args):
    refgene_fields = """
    bin
    name
    chrom
    strand
    txStart
    txEnd
    cdsStart
    cdsEnd
    exonCount
    exonStarts
    exonEnds
    score
    name2
    cdsStartStat
    cdsEndStat
    exonFrames
    """.split()
    #Skip the header lines, read in only the columns we need because some unnecessary columns can be millions of characters long
    reader = csv.DictReader(filter(lambda row: row[0]!='#',open(args.refgene,'r')), delimiter='\t', fieldnames=refgene_fields)
    out=[x for x in reader if x['chrom'] in chromosomes]
    sorted_out = natsorted(out, key=itemgetter('chrom'))
    headers = ['chrom','txStart','txEnd','name2','name','strand','exonCount','exonStarts','exonEnds']
    writer = csv.DictWriter(open(args.outfile,'w'), extrasaction='ignore',fieldnames=headers, delimiter='\t')
    writer.writerows(sorted_out)
