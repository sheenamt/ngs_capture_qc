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
    git clone git@github.com:sheenams/ngs_capture_qc.git
    cd ngs_capture_qc


execution
=========

The ``cap_qc.py`` script provides the different scripts used to process
the original probe file. To run ``cap_qc.py``::

    % ./cap_qc.py -h

The expected order of operations:
1. ./cap_qc.py filter_refseq 
   - this will validated preferred transcripts and create the filtered refseq file for use in CNV calling
2. ./cap_qc.py create_files
   - creates the following files:
     - clean bed (probes merged, deduplicated and annotated)
     - picard bed (probes in format required by Picard)
     - overall_summary (unique bases targeted, coding bases targeted, refseqs with at least 1 base targeted, probes outside of coding)
     - per refseq summary (total_bases_targeted,length_of_gene,fraction_of_gene_covered,exons_with_coverage)

Commands are constructed as follows. Every command starts with the
name of the script, followed by an "action" followed by a series of
required or optional "arguments". The name of the script, the action,
and options and their arguments are entered on the command line
separated by spaces. Help text is available for both the ``cap_qc``
script and individual actions using the ``-h`` or ``--help`` options::

    % cap_qc.py -h
    usage: cap_qc.py [-h] [-V] [-v] [-q]
    {help,create_files,filter_refseq,summarize_assay} ...

    Utilities for the ngs_capture_qc scripts

    positional arguments:
    {help,create_files,filter_refseq,summarize_assay}
    help                Detailed help for actions using `help <action>`
    create_files        Script to create specifically formatted files from the
                        original probes file
    filter_refseq       Filter a file containing the refGene annotation table,
                        limiting to
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

    % ./cap_qc.py -V
    0309.004ecac

unit tests
==========

Unit tests are implemented using the ``unittest`` module in the Python
standard library. The ``tests`` subdirectory is itself a Python
package that imports the local version (ie, the version in the project
directory, not the version installed to the system) of the ``munge``
package. All unit tests can be run like this::

     % ./testall
    ........................
    ----------------------------------------------------------------------
    Ran 7 tests in 0.155s

    OK

A single unit test can be run by referring to a specific module,
class, or method within the ``tests`` package using dot notation::

    % ./testone tests.test_utils
    .
    ----------------------------------------------------------------------
    Ran 1 test in 0.004s

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
