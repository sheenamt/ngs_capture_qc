"""
Create xlsx workbook from all output files

usage:

 capqc xlsmaker /path/to/summary/files

"""
import csv
import os

from xlsxwriter import Workbook

def build_parser(parser):
    parser.add_argument(
        'infiles', action='append', nargs='+',
        help='Input files')
    parser.add_argument(
        '-o', '--outfile',
        help='Output file', required=True)

def write_workbook(sheet_name, book, fname):
    """
    Write analysis file as sheet in workbook
    """
    sheet = book.add_worksheet(sheet_name)    
    Reader = csv.reader(open(fname, 'rU'), delimiter='\t')
    for rowx, row in enumerate(Reader):
        for colx, value in enumerate(row):
            sheet.write(rowx, colx, value)

def action(args):

    book = Workbook()
    (infiles, ) = args.infiles
    for fname in infiles:
        (f_path, f_name) = os.path.split(fname)
        (f_short_name, f_extension) = os.path.splitext(f_name)
        print(fname)
        write_workbook(f_short_name, book, fname)
    book.filename=args.outfile

    book.close()
