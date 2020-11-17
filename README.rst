==================================================================================
ngs_capture_qc: Python package for QC analysis of NGS capture probe/manifest files
==================================================================================

This project provides scripts for analyzing the probe files for an NGS capture
library and creating input files for various programs (picard, bedtools, etc)

.. contents:: Table of Contents

dependencies
============

* Python 3.3+
* Tested on Linux and OS X.
* bedtools >= 2.26

installation
============

Clone the project from the git repository::

    cd ~/src
    git clone git@gitlab.labmed.uw.edu:sheenams/ngs_capture_qc.git
    cd ngs_capture_qc
    python -m venv capqc-env
    source capqc-env/bin/activate
    pip install -r requirements.txt
    
execution
=========

The ``capqc`` script provides the different scripts used to process
the original probe file. To run ``capqc``::

    % ./capqc -h

The expected order of operations:
1. ./capqc refgene_to_bed [-h] refgene output
   - created BED format of refgene, if one isn't available already

create_files and summarize_assay expect refgene in bed format. 

2. ./capqc summarize_assay [-h] [--outdir OUTDIR] bed genes refgene_bed bedtools
    - overall_summary (unique bases targeted, coding bases targeted, refgenes with at least 1 base targeted, probes outside of coding)
    - preferred refgene summary (total_bases_targeted,length_of_gene,fraction_of_gene_covered,exons_with_coverage)
    - other refgene summary (total_bases_targeted,length_of_gene,fraction_of_gene_covered,exons_with_coverage)

3. ./capqc xlsxmaker [-h] -o OUTFILE infiles [infiles ...]
   - per refgene summary
   - other genes that are covered
   - overall summary
     
4. ./capqc filter_refgene [-h] refgene genes outfile
   - the refgene file input for this is NOT bed format
   - this will validated preferred transcripts and create the filtered refgene file for use in CNV calling

5. ./capqc capqc create_files [-h] probefile refgene_bed bedtools outdir
   - creates the following files:
     - clean bed (probes merged, deduplicated and annotated)
     - picard bed (probes in format required by Picard)

Commands are constructed as follows. Every command starts with the
name of the script, followed by an "action" followed by a series of
required or optional "arguments". The name of the script, the action,
and options and their arguments are entered on the command line
separated by spaces. Help text is available for both the ``cap_qc``
script and individual actions using the ``-h`` or ``--help`` options::

    % capqc -h
    usage: capqc [-h] [-V] [-v] [-q]
    {help,create_files,filter_refgene,refgene_to_bed,summarize_assay}
    ...

    Utilities for the ngs_capture_qc scripts

    positional arguments:
    {help,create_files,filter_refgene,refgene_to_bed,summarize_assay}
    help                Detailed help for actions using `help <action>`
    create_files        Script to create specifically formatted files from the
                        original probes file
    filter_refgene      Filter a file containing the refgene annotation table,
                        limiting to
    refgene_to_bed      Convert UCSC refgene.txt files to BED format, for use
                        in summarize_assay script
    xlsxmaker           Create xlsx workbook from all output files
    summarize_assay     Given probe reference file, list of preferred
                        transcripts and refgene.bed,

    optional arguments:
    -h, --help            show this help message and exit
    -V, --version         Print the version number and exit
    -v, --verbose         Increase verbosity of screen output (eg, -v is
                        verbose, -vv more so)
    -q, --quiet           Suppress output

versions
========

We use abbrevited git sha hashes to identify the software version::

    % ./capqc -V
    0309.004ecac

unit tests
==========

Unit tests are implemented using the ``unittest`` module in the Python
standard library. The ``tests`` subdirectory is itself a Python
package that imports the local version (ie, the version in the project
directory, not the version installed to the system) of the ``cap_qc``
package. All unit tests can be run like this::

     % ./testall
     ...................
     ----------------------------------------------------------------------
     Ran 19 tests in 4.624s
     
     OK


license
=======

Copyright (c) 2019 Sheena Todhunter

Released under the MIT License:

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
