"""
Convert UCSC refgene.txt files to BED format, for use in summarize_assay script
"""
 
import os
import sys 
 
def build_parser(parser):
    parser.add_argument('refseq', help='UCSC table browser download')
    parser.add_argument('output', help='Output file', default = sys.stdout)

 
def get_int_list(l):
    return [int(i) for i in l.strip(',').split(',')]
 
def get_string_list(a):
    return ','.join([str(i) for i in a])
 
def bed_key(d):
    return([d['chrom'], int(d['chromStart']), int(d['chromEnd'])])

def action(args):
    genes = {}
    for line in open(args.refseq, 'r'):
        ls = line.strip().split('\t')
        if len(ls)<16:
            sys.exit("File expected to have the following columns: 'bin', 'name', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds', 'score', 'name2', 'cdsStartStat', 'cdsEndStat', 'exonFrames'")
        if 'bin' in ls[0]:
            continue
        starts =  get_int_list(ls[9])
        stops = get_int_list(ls[10])
        lengths = get_string_list([stop-start for stop,start in zip(stops, starts)])        
        relstarts = get_string_list(starts) #[start - int(ls[4]) for start in starts])
        relends = get_string_list(stops) #[start - int(ls[4]) for start in starts])
        # For format see http://genome.ucsc.edu/FAQ/FAQformat.html#format1
        features = ['chrom','chromStart','chromEnd','name', 'refseq','exonCount','exonSizes','exonStarts','exonEnds']
        gene_entry = dict([('chrom', ls[2].strip('chr')),
                           ('chromStart', ls[4]),
                           ('chromEnd', ls[5]),
                           ('name', ls[12]),
                           ('refseq', ls[1]),
                           ('exonCount', ls[8]),
                           ('exonSizes', lengths),
                           ('exonStarts', relstarts),
                           ('exonEnds', relends)])
        refseq = ls[1]
        
        # Ensure that each refseq is only in the table once
        if (refseq not in genes 
            and 'NM' in refseq):
            genes[refseq] = gene_entry

    output=open(args.output,'w')
    for gene in sorted(genes.values(), key=bed_key):
        output.write('\t'.join([str(gene[f]) for f in features]) + '\n')
    output.close()

