"""
Parse the filtered refgene file and report chrm:start-stop for each gene
"""
 
import sys 
import csv
from operator import itemgetter

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
    sorted_out = sorted(out, key=itemgetter('name2'))
    headers = ['name2','name','chrom','txStart','txEnd']
    writer = csv.DictWriter(open(args.outfile,'w'), extrasaction='ignore',fieldnames=headers, delimiter='\t',lineterminator='\n')
    writer.writerows(sorted_out)
