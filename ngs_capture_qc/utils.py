import os
import gzip
import logging
import shutil
import sys
import numpy as np
from collections import defaultdict
from intervaltree import Interval, IntervalTree

try:
    import bz2
except ImportError as err:
    def bz2_open(filename, mode, *args, **kwargs):
        sys.exit(err)
else:
    bz2_open = bz2.open if hasattr(bz2, 'open') else bz2.BZ2File

log = logging.getLogger(__name__)


def cast(val):
    """Attempt to coerce `val` into a numeric type, or a string stripped
    of whitespace.

    """

    for func in [int, float, lambda x: x.strip(), lambda x: x]:
        try:
            return func(val)
        except ValueError:
            pass


def mkdir(dirpath, clobber=False):
    """
    Create a (potentially existing) directory without errors. Raise
    OSError if directory can't be created. If clobber is True, remove
    dirpath if it exists.
    """

    if clobber:
        shutil.rmtree(dirpath, ignore_errors=True)

    try:
        os.mkdir(dirpath)
    except OSError:
        pass

    if not os.path.exists(dirpath):
        raise OSError('Failed to create %s' % dirpath)

    return dirpath


class Opener(object):
    """Factory for creating file objects. Transparenty opens compressed
    files for reading or writing based on suffix (.gz and .bz2 only).

    Example::

        with Opener()('in.txt') as infile, Opener('w')('out.gz') as outfile:
            outfile.write(infile.read())
    """

    def __init__(self, mode='r', *args, **kwargs):
        self.mode = mode
        self.args = args
        self.kwargs = kwargs
        self.writable = 'w' in self.mode

    def __call__(self, obj):
        if obj is sys.stdout or obj is sys.stdin:
            return obj
        elif obj == '-':
            return sys.stdout if self.writable else sys.stdin
        else:
            openers = {'bz2': bz2_open, 'gz': gzip.open}
            __, suffix = obj.rsplit('.', 1)
            # in python3, both bz2 and gz libraries default to binary input and output
            mode = self.mode
            if sys.version_info.major == 3 and suffix in openers \
               and mode in {'w', 'r'}:
                mode += 't'
            opener = openers.get(suffix, open)
            return opener(obj, mode=mode, *self.args, **self.kwargs)

def check_probe_format(probes):
    """Check that the probes are in chrm|start|stop|annotation|strand format.
    Remove 'chr' if present"""
    assert len(probes.columns)>=5, "Five columns expected. Please format input file as chrm|start|stop|annotation|strand, without a header"
    #assert that chrm is in chromosome dictionary (ie, there is no header)
    if probes.iloc[0][0] not in chromosomes.keys():
        raise ValueError("Column 1 is not an obvious chromosome. Please format input file as chrm|start|stop|annotation|strand, without a header")
    elif not isinstance(probes.iloc[0][3], str):
        raise ValueError("Column 4 is not an obvious annotation. Please format input file as chrm|start|stop|annotation|strand, without a header")
    elif probes.iloc[0][4] not in ['-','+']:
        raise ValueError("Column 5 is not an obvious strand (-,+). Please format input file as chrm|start|stop|annotation|strand, without a header")
    elif isinstance(probes.iloc[0][1],str) or isinstance(probes.iloc[0][2],str) :
        raise ValueError("Column 2 and/or 3 is not an obvious start|stop position. Please format input file as chrm|start|stop|annotation|strand, without a header")

    #Drop all other columns
    probes=probes.iloc[:,:5]
    probes.columns=['chrom','start','stop','annotation','strand']
    probes=probes.astype({"chrom":str})
    #Drop chr if present
    probes['chrom']=probes['chrom'].str.replace('chr','')
    return probes

# Various files and data strctures specify chromosomes as strings
# encoding ints, like ('1', '2', ..., 'X'), sometimes as ints (1, 2,
# ... 'X'), and sometimes with a prefix ('chr1', 'chr2', ...,
# 'chrX'). `chromosomes` maps all three to the numeric representation.
chrnums = list(range(1, 23)) + ['X', 'Y']
chromosomes = {'chr{}'.format(c): c for c in chrnums}
chromosomes.update({str(c): c for c in chrnums})
chromosomes.update({c: c for c in chrnums})

class UCSCTable(object):
    '''A container class for the parsing functions, used in GenomeIntervalTree.from_table``.'''
    REF_GENE_FIELDS = ['bin', 'name', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds', 'score', 'name2', 'cdsStartStat', 'cdsEndStat', 'exonFrames']
    @staticmethod
    def REF_GENE(line):
        return dict(zip(UCSCTable.REF_GENE_FIELDS, line.split(b'\t')))

class IntervalMakers(object):
    '''A container class for interval-making functions, used in GenomeIntervalTree.from_table and GenomeIntervalTree.from_bed.'''


    @staticmethod
    def TX(d):
        return [Interval(int(d['txStart']), int(d['txEnd']), d)]

    @staticmethod
    def CDS(d):
        return [Interval(int(d['cdsStart']), int(d['cdsEnd']), d)]

    @staticmethod
    def EXONS(d):
        exStarts = d['exonStarts'].split(b',')
        exEnds = d['exonEnds'].split(b',')
        intron_count=int(d['exonCount'])-1
        exon_count=int(d['exonCount'])
        strand = d['strand']
        for i in range(exon_count):
            exon_d = d.copy()
            if strand == '+':
                exon_d['exonNum']=str(i+1)
            elif strand == '-':
                exon_d['exonNum']=str(exon_count-i)
            #Since interval trees are not inclusive of upper limit, add one to the exon end boundary
            yield Interval(int(exStarts[i]), int(exEnds[i])+1, exon_d)

            #Setup the intron info
            if i < intron_count:
            #Since interval trees are not inclsive of upper limit, add one to the intron start boundary and not to the end boundary
                intron_start=int(exEnds[i])+1
                intron_end=int(exStarts[i+1])
                new_d=d.copy()
                if new_d['strand']=='-':
                    new_d['intronNum']=str(intron_count - i)
                elif new_d['strand']=='+':
                    new_d['intronNum']=str(i+1)
                yield Interval(intron_start, intron_end, new_d)

def _fix(interval):
    '''
    Helper function for ``GenomeIntervalTree.from_bed and ``.from_table``.

    Data tables may contain intervals with begin >= end. Such intervals lead to infinite recursions and
    other unpleasant behaviour, so something has to be done about them. We 'fix' them by simply setting end = begin+1.
    '''
    if interval.begin >= interval.end:
        log.info("Interval with reversed coordinates (begin >= end) detected when reading data. Interval was automatically fixed to point interval [begin, begin+1).")
        return Interval(interval.begin, interval.begin+1, interval.data)
    else:
        return interval

class GenomeIntervalTree(defaultdict):
    '''
    The data structure maintains a set of IntervalTrees, one for each chromosome.
    It is essentially a ``defaultdict(IntervalTree)`` with a couple of convenience methods
    for reading various data formats.
    '''
    def __init__(self):
        super(GenomeIntervalTree, self).__init__(IntervalTree)

    def addi(self, chrom, begin, end, data=None):
        self[chrom].addi(begin, end, data)

    def __len__(self):
        return sum([len(tree) for tree in self.values()])

    @staticmethod
    def from_table(fileobj=None, url='http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz',
                    parser=UCSCTable.REF_GENE, mode='tx', decompress=None):
        '''
        Index the rows of UCSC tables into a ``GenomeIntervalTree`` 

        The table can be either specified as a ``fileobj`` (in which case the data is read line by line),
        or via an ``url`` (the ``url`` may be to a ``txt`` or ``txt.gz`` file either online or locally).
        The type of the table is specified using the ``parser`` parameter. This is a function that takes a line
        of the file (with no line ending) and returns a dictionary, mapping field names to values. This dictionary will be assigned
        to the ``data`` field of each interval in the resulting tree.

        Finally, there are different ways genes can be mapped into intervals for the sake of indexing as an interval tree.
        One way is to represent each gene via its transcribed region (``txStart``..``txEnd``). Another is to represent using
        coding region (``cdsStart``..``cdsEnd``). Finally, the third possibility is to map each gene into several intervals,
        corresponding to its exons (``exonStarts``..``exonEnds``).

        The mode, in which genes are mapped to intervals is specified via the ``mode`` parameter. The value can be ``tx``, ``cds`` and
        ``exons``, corresponding to the three mentioned possibilities.

        The ``parser`` function must ensure that its output contains the field named ``chrom``, and also fields named ``txStart``/``txEnd`` if ``mode=='tx'``,
        fields ``cdsStart``/``cdsEnd`` if ``mode=='cds'``, and fields ``exonCount``/``exonStarts``/``exonEnds`` if ``mode=='exons'``.

        The ``decompress`` parameter specifies whether the provided file is gzip-compressed.
        This only applies to the situation when the url is given (no decompression is made if fileobj is provided in any case).
        If decompress is None, data is decompressed if the url ends with .gz, otherwise decompress = True forces decompression.

        '''
        #Read in data from URL if file not provided
        if fileobj is None:
            data = urlopen(url).read()
            if (decompress is None and url.endswith('.gz')) or decompress:
                data = zlib.decompress(data, 16+zlib.MAX_WBITS)
            fileobj = BytesIO(data)

        interval_lists = defaultdict(list)

        #Setup the interval type
        if mode == 'tx':
            interval_maker = IntervalMakers.TX
        elif mode == 'cds':
            interval_maker = IntervalMakers.CDS
        elif mode == 'exons':
            interval_maker = IntervalMakers.EXONS
        elif getattr(mode, __call__, None) is None:
            raise Exception("Parameter `mode` may only be 'tx', 'cds', 'exons' or a callable")
        else:
            interval_maker = mode

        #Parse the genome data
        for ln in fileobj:
            if not isinstance(ln, bytes):
               ln = ln.encode()
            if ln.startswith(b'#'): #python2 needs:or ln.startswith('#'):
                continue
            ln = ln.strip()
            d = parser(ln)
            for interval in interval_maker(d):
                interval_lists[d['chrom']].append(_fix(interval))
                
        # Now convert interval lists into trees
        gtree = GenomeIntervalTree()
        for chrom, lst in getattr(interval_lists, 'iteritems', interval_lists.items)():
            gtree[chrom] = IntervalTree(lst)
        return gtree
        
    def __reduce__(self):
        t = defaultdict.__reduce__(self)
        return (t[0], ()) + t[2:]

