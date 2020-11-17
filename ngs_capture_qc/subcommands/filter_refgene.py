"""Filter a file containing the refgene annotation table, limiting to
preferred transcripts.

Overlapping genes will result in an error.
"""

import sys
import csv
import pprint
from itertools import groupby
from operator import itemgetter
from collections import namedtuple
import logging
import pandas as pd

from ngs_capture_qc.utils import chromosomes

log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('refgene',
                        help='RefSeq table broswer file, NOT BED FORMAT')
    parser.add_argument('genes', 
                        help='File defining preferred transcripts')
    parser.add_argument('outfile',
                        help='output file')

refgene_fields = """
bin
refgene
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
name
cdsStartStat
cdsEndStat
exonFrames
""".split()

def read_refgene(file):
    """Read open file-like object `file` containing annotations and
    return a sequence of dicts.

    The `annotations` file is available from the UCSC Genome Browser
    website:
    http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz

    To view the schema, go to https://genome.ucsc.edu/cgi-bin/hgTables
    --> choose track "refSeq Genes" --> choose table "refGene" -->
    "describe table schema"
    
    Skip the header (which starts with #bin) if present
    """
    
    return csv.DictReader(filter(lambda row: row[0]!='#',file), fieldnames=refgene_fields, delimiter='\t')


def check_overlapping(features):
    """Check for elements of `features` with overlapping ranges. In the
    case of overlap, print an informative error message and return
    names and positions of overlapping features.

    """
    features = features[:]
    overlapping = []
    for i in range(len(features)-1):
        prev_name, prev_start, prev_end = features[i]
        name, start, end = features[i+1]
        if prev_end >= start:
            overlap = ((prev_name, prev_start, prev_end), (name, start, end))
            overlapping.append(overlap)
            raise ValueError('overlapping features: ' + pprint.pformat(overlap))
    return overlapping

Node = namedtuple('Node', 'name start end left right')

def action(args):
    #Parse preferred transcripts, write coverage info
    transcripts=pd.read_csv(open(args.genes, 'r'), delimiter='\t')
    transcripts.columns=['Gene','RefSeq']
    transcripts['RefSeq']=transcripts['RefSeq'].apply(lambda x: str(x).split('.')[0])

    # read and filter the refgene file
    refgenes = read_refgene(open(args.refgene, 'r'))
    fieldnames = refgenes.fieldnames

    refgenes = [r for r in refgenes
                if r['chrom'] in chromosomes and r['name'] in transcripts['Gene'].values]

    # sort by chromosome, transcription start
    refgenes.sort(key=lambda row: (str(chromosomes[row['chrom']]), str(row['txStart'])))

    # group by gene and choose one transcript for each
    filtered_output = []
    for gene, grp in groupby(refgenes, itemgetter('name')):
        grp = list(grp)
        preferred = transcripts.loc[transcripts['Gene']==gene, 'RefSeq'].item()
        if preferred:
            keep = [r for r in grp if r['refgene'].split('.')[0] in preferred]
            if not keep:
                log.error('Error: %s has a preferred transcript of %s but only %s was found' %
                          (gene, preferred, ','.join(r['refgene'] for r in grp)))
                sys.exit(1)
            elif len(keep) > 1:
                log.warning('{} has more than one preferred transcript; using {}'.format(
                    gene, keep[0]['refgene']))
        else:
            log.warning('no preferred transcript for {}'.format(gene))
            keep = grp[:1]

        filtered_output.append(keep[0])

#    all transcripts are found among preferred transcripts
    # for k in filtered_output:
    #     print(k)
    #     print(transcripts[transcripts['RefSeq'].astype(str).str.contains(k['refgene'])])
    assert [transcripts[transcripts['RefSeq'].astype(str).str.contains(k['refgene'])] for k in filtered_output]

    # Check for overlapping genes and exit with an error if any are
    # found.
    overlapping = []
    rows = sorted(filtered_output, key=itemgetter('chrom'))
    for chrom, grp in groupby(rows, key=itemgetter('chrom')):
        grp = sorted(
            [(row['name'], int(row['txStart']), int(row['txEnd'])) for row in grp],
            key=itemgetter(1, 2)
        )
        overlapping.extend(check_overlapping(grp))

    if overlapping:
        log.error('Error: overlapping genes were found')
        sys.exit(1)


    writer = csv.DictWriter(open(args.outfile,'w'), fieldnames=fieldnames, delimiter='\t')
    writer.writerows(filtered_output)
