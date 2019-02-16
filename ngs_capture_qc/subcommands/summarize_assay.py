"""
Given probe reference file, list of preferred transcripts and refgene.bed, 
compute per-refgene and summary statistics and output any genes not covered as expected
Requires refgene data in bed format, available with UW genomes or use refgene_to_bed script to create 
"""
 
import sys
import subprocess
import csv
import os
import logging 
import pandas as pd
import numpy as np

from collections import defaultdict
from ngs_capture_qc.utils import chromosomes
if sys.version_info[0] < 3: 
    from StringIO import StringIO
else:
    from io import StringIO
log = logging.getLogger(__name__)

def build_parser(parser):
    parser.add_argument('--bed', required=True, help="Assay Reference bed file, sorted, with ^M removed from end of lines")
    parser.add_argument('--genes', required=True, help="Gene, RefSeq for assay")
    parser.add_argument('--refgene_bed', help="UCSC Refgene data in bed format")
    parser.add_argument('--bedtools', default='',help='Path to bedtools, accepts binary or singularity image')
    parser.add_argument('--outdir', required=False, help="Output directory for summary scripts")

class exonTracker:
    """
    Keeps track of a gene's exons.  When an interval is inserted, any relevant exons are covered.
    Initially, all exon intervals' values are covered=False
    """
    def __init__(self, exonStarts, exonEnds):
        self.exons = dict(((start,end),False) for start,end in zip(exonStarts,exonEnds))
        assert(len(exonStarts) == len(exonEnds))
    def insert(self, start, end):
        #Probes are merged, and therefore may cover multiple exons
        #probe: 1-100, exons: 1-20, 30-49, 59-99
        for exonStart, exonEnd in self.exons.keys():
            if (int(start) >= int(exonStart) and int(end) < int(exonEnd)) or (int(end) > int(exonStart) and int(end) <= int(exonEnd)) or (int(exonStart) > int(start) and int(exonEnd) < int(end)):
                self.exons[(exonStart,exonEnd)] = True
                
def calculate_total_covered(probes):
    '''calculate the total regions covered by using the merged probes file'''
    total_cov=0
    with open(probes, 'r') as p:
        for line in p:
            ls=line.split('\t')
            chrm,start,stop=ls[0:3]
            line_sum=int(stop)-int(start)
            total_cov += line_sum

    return total_cov


def action(args):

    if args.bedtools.endswith('img'):
        bedtools='singularity exec --bind {} --pwd {} {}'.format(os.getcwd(), os.getcwd(), args.bedtools)
    else:
        bedtools=args.bedtools

    out = args.outdir if args.outdir else ''
    refgenes = {}
    genes = {}

    refgene_header = ['chrom','chromStart','chromEnd','name', 'refgene','exonCount','exonSizes','exonStarts','exonEnds'] 
    probes_header = ['chrom', 'chromStart', 'chromEnd' ]
    genes_header = ['Gene', 'RefSeq']
    
    # 1) Read refGene.txt into the refgenes dictionary
    for line in csv.DictReader(open(args.refgene_bed, 'r'), delimiter='\t', fieldnames=refgene_header):
        refgene = line['refgene']
        name = line['name']
        # Dictionary-ize refgene.bed
        # Insert unseen refgenes into the dictionary; 
        # We asume that refgene only has ONE line per refgene
        exonStarts=[x for x in line['exonStarts'].split(',')]
        exonEnds=[x for x in line['exonEnds'].split(',')]
        if refgene not in refgenes:
            refgenes[refgene] = dict( [('name', name),
                                     ('refgene', line['refgene']),
                                     ('chrom', line['chrom'].strip('chr')),
                                     ('chromStart', int(line['chromStart'])),
                                     ('chromEnd', int(line['chromEnd'])),
                                     ('exonTracker', exonTracker(exonStarts, exonEnds)),
                                     ('bases_covered', 0)])
            
            #Sanity checks
            assert(len(exonStarts) == len(exonEnds))
            for start,end in zip(exonStarts, exonEnds):
                assert(int(start) < int(end))
                assert(int(line['chromStart']) <= int(start) and int(start) < int(line['chromEnd']))
                assert(int(line['chromStart']) < int(end) and int(end) <= int(line['chromEnd']))
        else:
            sys.stderr.write("RefSeq {} is listed twice in refGene!".format(line['refgene']))

    # 2) Using bedtools, calculate how many bases are actually covered for each gene
    intersection=os.path.join(out,'intersect_probes_refgene.txt')
    write_intersect=open(intersection, 'w')
    intersect_args = [x for x in bedtools.split(' ')]+['bedtools','intersect', '-wo' ,'-a', args.bed, '-b', args.refgene_bed]
    intersect = subprocess.call(intersect_args, stdout=write_intersect)
    write_intersect.close()

    # Parse that output, collecting the number of covered bases per-gene, and annotate refgenes dictionary
    # Note: Communicate returns (stdoutdata, stderrdata), stdout is a giant string, not an iterable
    # Also, the last line is just a newline, which must be skipped
    for line in open(intersection):
        ls = line.split('\t')
        #Find the NM_ column, can be different depending on the input file
        indices = [i for i, s in enumerate(ls) if 'NM_' in s.upper()]
        if len(indices)>1:
            sys.stderr.write("Refseq {} is listed twice in refGene!".format(line['refgene']))
        refgene = ls[indices[0]]          # We pick out the refgene of the gene from refGene that was matched
        overlap = int(ls[-1]) # The '-wo' switch from intersect_args put the amount of overlap here
        refgenes[refgene]['bases_covered'] += overlap
        refgenes[refgene]['exonTracker'].insert(int(ls[1]), int(ls[2]))


    # 4) Print per-refgene summary
    per_refgene_header = ['gene','refgene','total_bases_targeted','length_of_gene','fraction_of_gene_covered','exons_with_any_coverage','total_exons_in_gene']
    per_refgene_writer = csv.DictWriter(open(os.path.join(out, "per_refgene_summary.txt"), 'w'), fieldnames=per_refgene_header,  delimiter='\t', extrasaction='ignore')
    per_refgene_writer.writeheader()
    # While we're looping through refgenes, count the total bases,and refgenes covered
    total_coding_bases = 0
    gene_count = 0

    for gene in csv.DictReader(open(args.genes, 'r'), delimiter='\t', fieldnames=genes_header):
        transcript = gene['RefSeq'].split('.')[0]
        if transcript.upper()=='REFSEQ':
            continue
        try:
            #assert that the NM provided in the Pref_Trans is for the Gene they requested
            if refgenes[transcript]['name']!=gene['Gene']:
                outfields=wrong_ref(gene)
                outfields = dict([('gene', gene['Gene']), 
                                  ('refgene', gene['RefSeq']),
                                  ('total_bases_targeted', 'Incorrect RefSeq for this Gene'),
                                  ('length_of_gene','NA'),
                                  ('fraction_of_gene_covered','NA'),
                                  ('exons_with_any_coverage','NA'),
                                  ('total_exons_in_gene','NA')])
            else:
                gene['bases_covered']=refgenes[transcript]['bases_covered']
                #Only count this as a covered gene if it has coverage
                if gene['bases_covered'] > 0:
                    gene_count +=1
                exons = [exon for exon in refgenes[transcript]['exonTracker'].exons.values()].count(True)

                outfields = dict([('gene', gene['Gene']), 
                                  ('refgene', transcript),
                                  ('total_bases_targeted', gene['bases_covered']),
                                  ('length_of_gene',refgenes[transcript]['chromEnd'] - refgenes[transcript]['chromStart']),
                                  ('fraction_of_gene_covered',round(float(gene['bases_covered']) /
                                                                    float(refgenes[transcript]['chromEnd'] - refgenes[transcript]['chromStart']),3)),
                                  ('exons_with_any_coverage',exons),
                                  ('total_exons_in_gene',len(refgenes[transcript]['exonTracker'].exons))])
                total_coding_bases += gene['bases_covered']

        #If this refgene isn't found, we should state that, cleanly 
        except KeyError:
            outfields = dict([('gene', gene['Gene']), 
                              ('refgene', gene['RefSeq']),
                              ('total_bases_targeted', 'RefSeq not found'),
                              ('length_of_gene','NA'),
                              ('fraction_of_gene_covered','NA'),
                              ('exons_with_any_coverage','NA'),
                              ('total_exons_in_gene','NA')])
        genes[gene['Gene']] = outfields

    
    #Process all data in the refgene file that overlaps our probes
    for transcript,data in refgenes.items():
        #add data for probes we didn't plan to cover
        if data['name'] in genes.keys() :
            continue
        else:
            if data['bases_covered'] > 0:
                gene_count +=1
                exons = [exon for exon in data['exonTracker'].exons.values()].count(True)
                outfields = dict([('gene', data['name']), 
                                  ('refgene', transcript),
                                  ('total_bases_targeted', data['bases_covered']),
                                  ('length_of_gene',data['chromEnd'] - data['chromStart']),
                                  ('fraction_of_gene_covered',round(float(data['bases_covered']) /
                                                                    float(data['chromEnd'] - data['chromStart']),3)),
                                  ('exons_with_any_coverage',exons),
                                  ('total_exons_in_gene',len(data['exonTracker'].exons))])
                total_coding_bases += data['bases_covered']
                genes[data['name']] = outfields

    for gene,data in sorted(genes.items()):
        per_refgene_writer.writerow(data)



    #5)  Calculate total regions covered 
    total_bases = calculate_total_covered(args.bed)
    non_intersection=os.path.join(out,'non_intersect_probes_refgene.txt')
    write_non_intersect=open(non_intersection, 'w')
    non_intersect_args = [x for x in bedtools.split(' ')]+['bedtools','intersect', '-v' ,'-a', args.bed, '-b', args.refgene_bed]
    non_intersect = subprocess.call(non_intersect_args, stdout=write_non_intersect)
    write_non_intersect.close()

    # 6) Print overall summary
    overall = open(os.path.join(out, "overall_summary.txt"),'w')

    # Note: The total bases and exon counts are probably slightly overestimated, since refgenes can
    # overlap and share bases.  The number of overlapping bases and exons, however, are neglible
    # and cumbersome to calculate
    overall.write("{} unique bases were targeted\n".format(total_bases))
    overall.write("{} unique bases within gene boundaries were targeted\n".format(total_coding_bases))
    overall.write("{} unique refgenes had at least one base targeted\n".format(gene_count))

    data=open(non_intersection)
    if data:
        overall.write("The following probes did not intersect with transcription region of any UCSC gene:\n")
        for line in data:
            overall.write(line)


