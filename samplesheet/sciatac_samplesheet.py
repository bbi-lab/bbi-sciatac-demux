#!/usr/bin/env python3

#
# Notes:
#   o  the master file is bbi-sciatac-demux/samplesheet/sciatac_samplesheet.py

"""
Program: sciatac_samplesheet.py
Summary:
  This program reads a (front-end) CSV spreadsheet file, performs a variety of
  checks, 'fixes' sample names, and writes a 'back-end' samplesheet file for
  use with the processing pipeline. There are two back-end samplesheet formats:
  Andrew Hill's sci-ATAC pipeline samplesheet format and a JSON format file
  required by the BBI sci-ATAC pipeline.

Input (front-end) samplesheet format:
  o  the input samplesheet file is a CSV format spreadsheet file (use a
     spreadsheet program to create the CSV file)
  o  checks rows
       o  checks the first eight cells in each row and trims off additional cells
       o  ignores rows with all empty cells
       o  reports an error if there are both empty and non-empty cells amongst
          the first six (the last two cells in a non-header row may be empty)
  o  the first row is a header with the following required columns
       o  n5 barcode identifier with possible values: 'n5_wells' or 'n5_indexes'
       o  n7 barcode identifier with possible values: 'n7_wells' or 'n7_indexes'
       o  p5 barcode identifier with possible values: 'p5_wells' or 'p5_indexes'
          or 'p5_columns'. If the p5 identifier is p5_columns, then the p7
          identifier must be p7_rows.
       o  p7 barcode identifier with possible values: 'p7_wells' or 'p7_indexes'
          or 'p7_rows'. If the p7 identifier is p7_rows, then the p5
          identifier must be p5_columns.
       o  sample name identifier with value: 'sample_name'
       o  genome label identifier with value: 'genome'
       o  peak_group (called peaks are merged by peak group)
       o  peak_file (a bed file of peaks; will be merged with called peaks if the
          peak_group value has non-zero length)
       o  external_sample_name the name of the sample supplied by the sample
          submitter
     and the following optional columns
       o  tissue the name of the tissue from which the sample was collected
       o  wrap_group the name of the group to which the results will
          be distributed. This information is used by the bbi-sciatac-wrap
          script.
  o  the column order is arbitrary (but must be consistent within the file)
  o  the header values do not depend on case
  o  ranges: index, well, column, and row range values are separated by either
     '-' or ':', your choice (mixes are allowed)
     examples:
       o  3-5
       o  3:5
       o  3-5,7:12
  o  multiple barcode values: multiple indexes, wells, columns, and rows (and
     ranges) are separated by ','
     examples:
       o  indexes: 1-5,7,89-96
       o  wells:   P1-A10,P1-A11
       o  rows:    E,F,G
       o  columns: 3,7,12
  o  sample names:
       o  begin with an alphanumeric characters: a-z, A-Z, and 0-9
       o  allowed characters are alphabetic (a-z and A-Z), numeric (0-9), and
          period '.'
       o  this program converts other characters in the sample name to '.',
          and then checks for sample name degeneracy. If there is, the program
          exits immediately.
  o  external sample names
       o  cannot contain the following characters: backspace, form feed,
          newline, carriage return, double quote, or backslash. Contiguous
          tabs are converted to a single space.
  o  tissue
       o  cannot contain the following characters: backspace, form feed, 
          newline, carriage return, double quote, or backslash. Contiguous
          tabs are converted to a single space.
  o  genomes:
       o  recognized genome names are listed in the variable 'genome_name_list'
          in this program's code. If a samplesheet genome name is not in the
          list, sciatac_samplesheet.py gives a warning in case the name is
          mis-spelled.
       o  genome is the name of the organism that was sequenced. This
          identifies the files required to analyze the reads. The genome string
          is passed to the processing pipeline, and the pipeline uses it to
          find required files. Available genomes are defined by the 'name'
           field in the bbi-sciatac-analyze/genomes.json file.
  o  peak groups
       o  the called peaks of samples with the same peak group name are merged
          prior to downstream analysis
       o  a sample peak group string (cell) may be empty
       o  peak group names consist of alphabetic, positive integers, and
          underscore characters
       o  each sample must have a peak_group or a peak_file or both
       o  if both a peak group and peak file are specified, the peaks in
          the peak bed file are merged with the called peaks in the group
          during the pipeline run
  o  peak files
       o  an absolute path to a peak bed file (it must begin with '/')
       o  a sample peak file string (cell) is empty when no peak
          file is wanted for the sample
       o  each sample must have a peak_group or a peak_file or both
       o  if both a peak group and peak file are specified, the peaks in
          the peak bed file are merged with the called peaks in the group
          during the pipeline run
  o  wrap_group
       o  valid names consist of lower and upper case alphabetic
          characters, numerals, and the symbols '.', '_', and '-'. Cells may
          be empty in which case those samples are excluded from the wrap.
  o  wells:
       o  samplesheet wells are converted to indexes where indexes refer to
          physical wells, which are in the order used in Andrew's pipeline;
          that is, N7 and P7 indexes increase by column number along each row
          and N5 and P5 indexes increase by row letter down each column
       o  wells that include a plate identifier can have either '-' or ':'
          separating the plate and well, your choice (mixes within a range are
          allowed)
          examples:
            o  P1-A10
            o  P1:A10
       o  wells are given as single or ranges of wells separated by commas.
          examples:
            o  A01-A12,B01:B12
            o  P1:A01-P1:A12,P2-A01:P2-A12
       o  when the same set of wells is used from all four plates for a sample,
          use 'P*' for the plate number
          example:
            o  P*:A01-P*:A04 is the same as P1:A01-P1:A04,P2:A01-P2:A04,P3:A01-P3:A04,P4:A01-P4:A04
       o  wells names do not depend on case
          examples:
            o  P1-A10
            o  p1-a10
  o  indexes:
       o  indexes are given as single or ranges of integer well indexes
          separated by commas.
          example:
            o  3,5-8,89:96
       o  indexes refer to physical wells in the order used in Andrew's
          pipeline and range from 1 to 384. Indexes 1 to 96 refer to
          plate 1, 97 to 192 to plate 2, and so on, for P5, P7, N5, and
          N7.
  o  p7 rows:
       o  u-titer plate rows used for PCR reactions
       o  rows are given as single or ranges of rows separated by commas
          examples:
            o  E,F,G
            o  D-F,H
  o  p5 columns:
       o  u-titer plate columns used for PCR reactions
       o  columns are given as single or ranges of columns separated
          by commas
          examples:
            o  3,4,5
            o  3-5,10-12
  o  p7 rows and p5 columns:
       o  by default, the i-th p7 row is paired with the i-th p5 column.
          No other combinations of specified rows and columns are considered.
          Consequently, the order of the rows and columns must be correct.
          Technical notes:
            o  by default, in order to restrict the PCR index combinations to
               one P7 row and one P5 column, this program expands the number of
               spreadsheet rows per sample by the number of P7 rows/P5 columns
               pairs given for the sample in the input CSV file. For example, if
               for a sample the p7_rows specification is 'A-B' and the
               p5_columns is '6,5', this program generates (internally) two rows
               for the sample, the first with PCR row column pair A/6 and the
               second with PCR row/column pair B/5.
            o  if the --no_expand_pcr_rows_columns option is set, for each
               sample, the ATAC-seq pipeline takes all combinations of the P7
               and P5 indices specified for the sample. Consequently, if the
               P7 rows were given as A-B and the P5 columns were given as 1-2
               in the input CSV file, and this program passed the corresponding
               indices directly to the pipeline, the pipeline would accept PCR
               index pairs for row/column pairs A/1, A/2, B/1, and B/2.
  o  p7 wells and p5 wells
       o  this program converts wells to indices and passes the indices to
          the pipeline. The pipeline uses combinations of all p7 and p5
          wells. The order of the wells is unimportant.
  o  name the barnyard sample 'Barnyard' for compatibility with the
     experiment dashboard

  Example samplesheet file:

  N7_indexes,N5_wells,P7_rows,P5_columns,sample_name,genome,peak_group,peak_file,external_sample_name,tissue,wrap_group
  1:96,P1-A01:P1-H01,"E,F,G","1,2,3",sample.1,mouse,group_1,,exsample.1,,Smith
  1:96,P1-A02:P1-H02,"E,F,G","1,2,3",sample.2,human,group_2,,exsample.2,,Jones
  1:96,P1-A03:P1-H03,"E,F,G","1,2,3",sample.3,mouse,group_1,,exsample.3,,Smith
  1:96,P1-A04:P1-H04,"E,F,G","1,2,3",sample.4,human,group_2,/home/me/my_peaks.bed,exsample.4,,Jones
  1:96,P1-A05:P1-H05,"E,F,G","1,2,3",sample.5,human,group_2,,exsample.5,,Jones
  1:96,P1-A06:P1-H06,"E,F,G","1,2,3",sample.6,mouse,group_1,,exsample.6,,Smith
  1:96,P1-A07:P1-H07,"E,F,G","1,2,3",sample.7,mouse,group_1,,exsample.7,,Smith
  1:96,P1-A08:P1-H08,"E,F,G","1,2,3",sample.8,mouse,group_1,,exsample.8,,Smith
  1:96,P1-A09:P1-H09,"E,F,G","1,2,3",sample.9,mouse,group_1,,exsample.9,,Smith
  1:96,P1-A10:P1-A10,"E,F,G","1,2,3",sample.10,human,group_2,,exsample.10,,Jones
  1:96,P1-A11:P1-H11,"E,F,G","1,2,3",sample.11,human,group_2,,exsample.11,,Jones
  1:96,P1-A12:P1-H12,"E,F,G","1,2,3",sample.12,barnyard,group_3,,exsample.12,,
  
  Notes:
    o  the p7_row and p5_column value sets are enclosed in quotes
       in this example because the sets have commas.  Do not use
       quotes in the spreadsheet program cells; the spreadsheet
       program adds them when it writes the CSV file.

Command line options:

  Option                   Description
  ------                   -----------
  --use_all_barcodes       By default, the BBI sci-ATAC processing pipeline uses only
                           informative barcodes to assign reads to samples. In this
                           description, a barcode is one of the four sequences N5, N7,
                           P5, P7 (for ligation and PCR reactions). These four
                           barcodes identify the read. If a barcode uses the same
                           index set for all samples, it is considered to be
                           uninformative for the purpose of identifying the sample,
                           so it is not used to identify the sample. In the example
                           samplesheet above, the N7, P5, and P7 index sets are the
                           same for all samples so only the N5 indexes are used to
                           identify samples. This is safe unless the sequencing run
                           includes samples that are not present in the samplesheet
                           and those additional samples are distinguished by barcode
                           index sets that appear to be uninformative based on the
                           samplesheet information. Use this option when there are
                           such additional samples in the sequencing run.

Notes:
  o  the command line arguments -r, -l, -w, -t, -s, --use_all_barcodes, and
     --hash_file are used only with the JSON samplesheet format.


For help with command line parameters, run

sciatac_samplesheet.py -h

which gives

usage: sciatac_samplesheet.py [-h] [-i INPUT] [-o OUTPUT] [-f {json,index}]
                              [-r RUN_DIR] [-l {2,3}] [-w {96,384}] [-t]
                              [-s {n5,n7}] [--use_all_barcodes] [-e] [-d] [-v]

A program to convert sci-ATAC CSV samplesheet to pipeline samplesheet.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input CSV samplesheet filename (required string).
  -o OUTPUT, --output OUTPUT
                        Output samplesheet filename (required string).
  -f {json,index}, --format {json,index}
                        Output file format (default: 'json') (optional
                        string).
  -r RUN_DIR, --run_dir RUN_DIR
                        Illumina run directory path (optional string).
  -l {2,3}, --level {2,3}
                        Two or three level sci-ATAC-seq experiment (default:
                        3) (optional integer).
  -w {96,384}, --number_wells {96,384}
                        Number of barcode set wells (default: 384) (optional
                        integer).
  -t, --tn5_barcodes    Tn5 has barcodes (optional flag).
  -s {n5,n7}, --sample_identifier {n5,n7}
                        Ligation barcode that identifies the sample (default:
                        'n5') used to check for duplicates (optional string).
  --hash_file {path}    Gives the full path to the sciPlex hash read index file
                        and enables hash read processing in the
                        bbi-sciatac-analyze pipeline.
  --use_all_barcodes    Use all barcodes to demultiplex fastq files. By
                        default, uninformative barcodes are not used to
                        identify samples when demultiplexing fastq files
                        (optional flag).
  -e, --template        Write template samplesheet file
                        ('samplesheet.template.csv') with standard column
                        formats and exit (optional flag).
  -d, --documentation   Display documentation and exit (optional flag).
  -v, --version         Give program and JSON output file versions and exit
                        (optional flag).
"""


import sys
import re
import csv
import json
import argparse

#
# Samplesheet JSON file version.
#
program_version = '4.2.0'
json_file_version = '3.2.0'

#
# List of recognizable genome names.
# This program issues a warning if a samplesheet genome is
# not in this list.
#
genome_name_list = [
  'arabidopsis',
  'barnyard',
  'bat',
  'cat',
  'chicken',
  'corn',
  'cow',
  'cynomolgus',
  'dog',
  'drosophila',
  'duck',
  'elephant',
  'horse',
  'human',
  'macaque',
  'mouse',
  'opossum',
  'pig',
  'rabbit',
  'rat',
  'snake',
  'worm',
  'zebrafish',
  'hg19',
  'mm19',
  'hg19_mm9'
]


#
# List of recognizable CSV column header names.
# These are used to check labels in the file.
# The n5, n7, p5, and p7 names are assigned to
# individual strings for error reporting.
#
n5_column_names = 'n5_wells n5_indexes'
n7_column_names = 'n7_wells n7_indexes'
p5_column_names = 'p5_wells p5_indexes p5_columns'
p7_column_names = 'p7_wells p7_indexes p7_rows'

column_header_name_list = [ 'sample_name', 'genome', 'peak_group', 'peak_file', 'external_sample_name', 'tissue', 'wrap_group' ]
column_header_name_list.extend( n5_column_names.split() )
column_header_name_list.extend( n7_column_names.split() )
column_header_name_list.extend( p5_column_names.split() )
column_header_name_list.extend( p7_column_names.split() )

columns_required_list = [ 'n5', 'n7', 'p5', 'p7', 'sample_name', 'genome', 'peak_group', 'peak_file', 'external_sample_name' ]
columns_optional_list = [ 'tissue', 'wrap_group' ]

#
# Columns that may have empty cells. This is used to prevent
# errors when testing for empty cells in check_rows().
#
column_allow_empty_cell = [ 'peak_group', 'peak_file', 'external_sample_name', 'tissue', 'wrap_group' ]

def display_documentation():
  print( __doc__ )
  return( 0 )


def clean_string(instring):
  # Remove leading and trailing whitespace.
  outstring = instring.strip()
  # Convert contiguous tabs to a space.
  outstring = re.sub( r'[\t]+', ' ', outstring )
  return( outstring )


def is_printable(string):
  # Check for non-printable characters.
  if(not string.isprintable()):
    return(False)
  return( True )


def has_space(string):
  if( re.search( r'[ ]', string ) ):
    return( True )
  return( False )


#
# Trim off leading and trailing whitespace and convert
# contiguous tabs to a space.
#
def clean_samplesheet_data( column_name_list, samplesheet_row_list ):
  errorFlag = False
  for irow, row_elements in enumerate(samplesheet_row_list):
    for icol in range( len( row_elements ) ):
      clean_string( row_elements[icol] )
      if( not is_printable( row_elements[icol] ) ):
        print('Error: non-printable character(s) in row %d column %d' % ( irow, icol + 1 ))
        errorFlag = True
  if( errorFlag ):
    sys.exit( -1 )


def check_args( args ):
  error_string = ''
  if( args.tn5_barcodes and ( args.number_wells == 384 or args.level == 3 ) ):
    error_string += '  tn5_barcodes requires level == 2 and number_wells == 96\n'
  if( len( error_string ) > 0 ):
    print( 'Command line argument errors:' )
    print( error_string, file=sys.stderr )
    sys.exit( -1 )
  return( 0 )


def parse_header_column_name( string_in, column_name_list, error_string ):
  """
  Split column header name into a 'type' and a 'format' and store as dictionary in column_name_list.
  """
  if( not string_in.lower() in column_header_name_list ):
    error_string += '  %s' % ( string_in )
    return( column_name_list, error_string )
  string_in = string_in.lower()
  if( string_in == 'sample_name' ):
    column_name_dict = { 'type': 'sample_name', 'format': None }
  elif( string_in == 'genome' ):
    column_name_dict = { 'type': 'genome', 'format': None }
  elif( string_in == 'peak_group' ):
    column_name_dict = { 'type': 'peak_group', 'format': None }
  elif( string_in == 'peak_file' ):
    column_name_dict = { 'type': 'peak_file', 'format': None }
  elif( string_in == 'external_sample_name' ):
    column_name_dict = { 'type': 'external_sample_name', 'format': None }
  elif( string_in == 'tissue' ):
    column_name_dict = { 'type': 'tissue', 'format': None }
  elif( string_in == 'wrap_group' ):
    column_name_dict = { 'type': 'wrap_group', 'format': None }
  else:
    mobj = re.match( r'([np][57])_(wells|indexes|rows|columns)', string_in )
    column_name_dict = { 'type': mobj.group( 1 ), 'format': mobj.group( 2 ) }
  column_name_list.append( column_name_dict )
  return( column_name_list, error_string )


def check_header_column_names( column_name_list ):
  """
  Check column header names.
  Notes:
    o  check that required columns occur
    o  check that each allowed column type occurs once
    o  check that if either p5 or p7 are specified by columns and rows, then both are.
  """
  columns_allowed = {}
  for column_name in columns_required_list:
    columns_allowed[column_name] = 0
  for column_name in columns_optional_list:
    columns_allowed[column_name] = 0
  for column_name_dict in column_name_list:
    columns_allowed[column_name_dict['type']] += 1
  errorFlag = 0
  # Check for n5, n7, p5, and p7 specification columns.
  for column_name in columns_allowed:
    if( columns_allowed[column_name] == 0 and column_name not in columns_optional_list ):
      print( 'Error: column for \'%s\' missing.' % ( column_name ), file=sys.stderr )
      if( column_name == 'n5' ):
        print( '  acceptable n5 header values: %s' % ( n5_column_names ), file=sys.stderr )
      elif( column_name == 'n7' ):
        print( '  acceptable n7 header values: %s' % ( n7_column_names ), file=sys.stderr )
      elif( column_name == 'p5' ):
        print( '  acceptable p5 header values: %s' % ( p5_column_names ), file=sys.stderr )
      elif( column_name == 'p7' ):
        print( '  acceptable p7 header values: %s' % ( p7_column_names ), file=sys.stderr )
      errorFlag = 1
    elif( columns_allowed[column_name] > 1 ):
      print( 'Error: column for \'%s\' occurs %d times.' % ( column_name, columns_allowed[column_name] ), file=sys.stderr )
      errorFlag = 1
  if( errorFlag ):
    sys.exit( -1 )
  p5_col = False
  p7_row = False
  for column_name_dict in column_name_list:
    if( column_name_dict['type'] == 'p5' and column_name_dict['format'] == 'columns' ):
      p5_col = True
    if( column_name_dict['type'] == 'p7' and column_name_dict['format'] == 'rows' ):
      p7_row = True
  if( p5_col != p7_row ):
    print( 'Error: p5 is %sin \'columns\' format but p7 is %sin \'rows\' format' % ( '' if p5_col else 'not ', '' if p7_row else 'not ' ), file=sys.stderr )
    sys.exit( -1 )
  return( 0 )


def parse_header( row_header ):
  """
  Convert column header (row) into a list of column name dictionaries.
  The dictionary has the elements
    key     value description
    type    entry type name: n5, n7, p5, p7, sample_name, genome, peak_group, peak_file, external_sample_name, tissue, wrap_group
    format  barcode format values: wells, indexes, rows, columns, None (see column_header_name_list for allowed combinations of type and format)
  """
  column_name_list = []
  error_string = ''
  for str in row_header:
    column_name_list, error_string = parse_header_column_name( str, column_name_list, error_string )
  if(len(error_string) > 0):
    print('Error: invalid header label(s): %s' % (error_string))
    sys.exit(-1)
  check_header_column_names( column_name_list )
  return( column_name_list )


def well_to_index( plate, row, column, across_row_first=True, element_coordinates=[None,None] ):
  """
  Convert a well specification to a plate index in the range P1:A01=1 to P4:H12=384.

  Args:
    plate              integer plate number between 1 and 4
    row                character row (A-H)
    column             integer column number (1-12)
    across_row_first   bool index increases by one as a row is traversed; that is,
                       moving from column to column along row
  Returns:
    index: an integer well index (1-384)

  """
  row = row.lower()
  if( plate < 1 or plate > 4  or
      not row in [ 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h' ] or
      column < 1 or column > 12 ):
    print( 'Error: spreadsheet cell: %s%s:  bad well values: plate: %d  row: %s  col: %d' % ( element_coordinates[0], element_coordinates[1], plate, row, column ), file=sys.stderr )
    sys.exit( -1 )
  irow = ord( row ) - ord( 'a' )
  icol = column - 1
  if( across_row_first ):
    well_index = irow * 12 + icol + 1
  else:
    well_index = icol * 8 + irow + 1
  return( well_index + ( plate - 1 ) * 96 )


def pad_well_col(well_col, zero_pad, id_length):
  if zero_pad:
      template = '%%0%sd' % id_length
  else:
      template = '%s'
  col_id = template % (well_col)
  return col_id


def index_to_well( well_index, across_row_first ):
  nrow = 8
  ncol = 12
  ipl = int( well_index / 96 )
  i96 = well_index - ipl * 96
  if across_row_first:
      well_row = chr(65 + int(i96 / ncol))
      well_col = (i96 % ncol) + 1
  else:
      well_row = chr(65 + (i96 % nrow))
      well_col = int(i96 / nrow) + 1

#    well_id = 'P%d-%s%s' % (ipl + 1, well_row, pad_well_col(well_col, zero_pad_col, id_length))
  well_id = '%s%s' % ( well_row, pad_well_col( well_col, True, 2 ) )

  return( well_id )


def check_index_list( index_list, element_coordinates = [ None, None ] ):
  """
  Check and clean index list
    o  check for duplicate indexes
    o  remove duplicate indexes
  """
  index_dict = {}
  for i in index_list:
    index_dict.setdefault( i, 0 )
    index_dict[i] += 1
  duplicate_list = []
  for i in index_dict.keys():
    if( index_dict[i] > 1 ):
      duplicate_list.append( str( i ) )
  if( len( duplicate_list ) > 0 ):
    print( 'Warning: spreadsheet cell: %s %s: duplicate well index(es): %s' % ( element_coordinates[0], element_coordinates[1], ' '.join( sorted(duplicate_list) ) ), file=sys.stderr )
    print( '         These indexes are duplicated within the reported spreadsheet cell, which may be intentional.')
  return( list( set( index_list ) ) )


def make_index_string( index_list ):
  """
  Convert a list of (integer) barcode well indexes to an index string where
    o  repeated indexes are dropped; that is, keep only distinct indexes
    o  sequences of counting numbers are expressed as ranges, for example, 5 6 7 8 9 => 5-9
    o  indexes and index ranges are separated by commas
  """
  index_string = ''
  index_list.sort()
  index_prev = None
  index1 = None
  for i in index_list:
    if( index_prev ):
      if( i == index_prev ):
        continue
      elif( i > index_prev + 1 ):
        if( len( index_string ) > 0 ):
          index_string += ','
        if( index_prev > index1 ):
          index_string += '%d-%d' % ( index1, index_prev )
        else:
          index_string += '%d' % ( index_prev )
        index1 = i
    else:
      index1 = i
    index_prev = i
  # last index in list
  if( len( index_string ) > 0 ):
    index_string += ','
  if( index_prev > index1 ):
    index_string += '%d-%d' % ( index1, index_prev )
  else:
    index_string += '%d' % ( index_prev )
  return( index_string )


def parse_indexes( string_in, max_index = 96, element_coordinates = [ None, None ] ):
  """
  Convert an index specification to an index string.
  Acceptable index specifications include
    o  single integer index
         o  76
    o  range of integer indexes, example
         o  5-9
    o  single and/or ranges of indexes separated by commas, examples
         o  3,6
         o  3,5-9
         o  2-4,9-18
         o  2:4,9:18
  """
  index_list = []
  string_in = re.sub( r'\s', '', string_in )
  for index_range in string_in.split( ',' ):
    mobj = re.match( r'([0-9]+)([-:]([0-9]+))?$', index_range )
    if( not mobj ):
      print( 'Error: spreadsheet cell: %s%s: bad index range \'%s\'' % ( element_coordinates[0], element_coordinates[1], index_range ), file=sys.stderr )
      sys.exit( -1 )
    index1 = int( mobj.group( 1 ) )
    if( index1 < 1 or index1 > max_index ):
      print( 'Error: spreadsheet cell: %s%s: bad index range \'%s\'' % ( element_coordinates[0], element_coordinates[1], index_range ), file=sys.stderr )
      sys.exit( -1 )
    index2 = int( mobj.group( 3 ) ) if mobj.group( 2 ) else index1
    if( index2 < 1 or index2 > max_index ):
      print( 'Error: spreadsheet cell: %s%s: bad index range \'%s\'' % ( element_coordinates[0], element_coordinates[1], index_range ), file=sys.stderr )
      sys.exit( -1 )
    if( index2 < index2 ):
      print( 'Error: spreadsheet cell: %s%s: bad index range \'%s\'' % ( element_coordinates[0], element_coordinates[1], index_range ), file=sys.stderr )
      sys.exit( -1 )
    for i in range( index1, index2 + 1 ):
      index_list.append( i )
  return( check_index_list( index_list, element_coordinates ) )


def parse_wells( string_in, across_row_first=True, max_index = 96, element_coordinates = [ None, None ] ):
  """
  Convert a well specification to an index string.
  Acceptable well specifications include
    o  a single well without a plate specified (implicit plate=1), examples
         o  A5
         o  A05
    o  a single well with a plate specified
         o  P1-A5
         o  P1:A5
    o  the same well from all four plates
         o  P*-A5
         o  P*:A5
    o  range of wells without plates specified,
       Note: the range of indices depends on whether the reaction type is n5/p5 or n7/p7.
         o  A9-B9
         o  A9:B9
    o  range of wells with plates specified,
         o  P1-A5:P1-A12
         o  P1:A5-P1:A12
         o  P1:A5:P1:A12
         o  P1-A5-P1-A12
    o  the same range of wells from all four plates
         o  P*-A5:P*-A10
    o  single and/or ranges of wells separated by commas
         o  A10,P1-B5:P1-B10,C7
  """
  index_list = []
  string_in = re.sub( r'\s', '', string_in )
  for well_range in string_in.split( ',' ):
    expand_plate_flag = False
#    mobj = re.match( r'(P([1-4])-)?([A-H])([01]?[0-9])([-:](P([1-4])-)?([A-H])([01]?[0-9]))?$', well_range )
    mobj = re.match( r'([pP][0]?([1-4*])[-:])?([a-hA-H])([0]?[1-9][0-2]?)([-:]([pP][0]?([1-4*])[-:])?([a-hA-H])([0]?[1-9][0-2]?))?$', well_range )
    if( not mobj ):
      print( 'Error: spreadsheet cell: %s%s: bad well or well range \'%s\'' % ( element_coordinates[0], element_coordinates[1], well_range ), file=sys.stderr )
      sys.exit( -1 )
    #
    # first well
    row1 = mobj.group( 3 )
    col1 = int( mobj.group( 4 ) )
    if( col1 < 1 or col1 > 12 ):
      print( 'Error: spreadsheet cell: %s%s: bad well: \'%s\'' % ( element_coordinates[0], element_coordinates[1], string_in ), file=sys.stderr )
      sys.exit( -1 )
    # is plate specified?
    if( mobj.group( 1 ) ):
      if( mobj.group( 2 ) == '*' ):
        plate1_list = [ 1, 2, 3, 4 ]
      else:
        plate1_list = [ int( mobj.group( 2 ) ) ]
    else:
      plate1_list = [ 1 ]
    #
    if( mobj.group( 5 ) ):
      if( ( mobj.group( 2 ) == None ) != ( mobj.group( 7 ) == None ) ):
        print( 'Error: spreadsheet cell: %s%s: either both or neither well in a range must have plates specified: \'%s\'' % ( element_coordinates[0], element_coordinates[1], string_in ), file=sys.stderr )
        sys.exit( -1 )
      if( ( mobj.group( 2 ) == '*' ) != ( mobj.group( 7 ) == '*' ) ):
        print( 'Error: spreadsheet cell: %s%s: either both or neither well in a range must have plates specified as \'*\': \'%s\'' % ( element_coordinates[0], element_coordinates[1], string_in ), file=sys.stderr )
        sys.exit( -1 )

      # second well, if this is a range
      row2 = mobj.group( 8 )
      col2 = int( mobj.group( 9 ) )
      if( col2 < 1 or col2 > 12 ):
        print( 'Error: spreadsheet cell: %s%s: bad well: \'%s\'' % ( element_coordinates[0], element_coordinates[1], string_in ), file=sys.stderr )
        sys.exit( -1 )
      # is plate specified?
      if( mobj.group( 6 ) ):
        if( mobj.group( 7 ) == '*' ):
          plate2_list = [ 1, 2, 3, 4 ]
        else:
          plate2_list = [ int( mobj.group( 7 ) ) ]
      else:
        plate2_list = [ 1 ]
      #
    else:
      plate2_list = plate1_list
      row2 = row1
      col2 = col1
    #
    for plate1, plate2 in zip( plate1_list, plate2_list ):
      index1 = well_to_index( plate1, row1, col1, across_row_first, element_coordinates )
      index2 = well_to_index( plate2, row2, col2, across_row_first, element_coordinates )
      if( index2 < index1 ):
        print( 'Error: spreadsheet cell: %s%s: bad well range: \'%s\'' % ( element_coordinates[0], element_coordinates[1], string_in ), file=sys.stderr )
        sys.exit( -1 )
      for i in range( index1, index2 + 1 ):
        index_list.append( i )
  return( check_index_list( index_list, element_coordinates ) )


def parse_rows( string_in, element_coordinates = [ None, None ] ):
  """
  Convert a row specification to an index string.
  Acceptable row specifications include
    o  single row
         o  B
    o  row range
         o  E-G
    o  single and/or ranges of rows separated by commas
         o  E-F,H
  """
  index_list = []
  string_in = re.sub( r'\s', '', string_in )
  for row_range in string_in.split( ',' ):
    mobj = re.match( r'([a-hA-H])([-:]([a-hA-H]))?$', row_range )
    if( not mobj ):
      print( 'Error: spreadsheet cell: %s%s: bad row or row range \'%s\'' % ( element_coordinates[0], element_coordinates[1], row_range ), file=sys.stderr )
      sys.exit( -1 )
    row1 = mobj.group( 1 )
    row1_index = well_to_index( 1, row1, 1, True, element_coordinates )
    row2_index = row1_index
    if( mobj.group( 2 ) ):
      row2 = mobj.group( 3 )
      row2_index = well_to_index( 1, row2, 1, True, element_coordinates )
      if( row2_index < row1_index ):
        print( 'Error: spreadsheet cell: %s%s: bad row range: \'%\'' % ( element_coordinates[0], element_coordinates[1], string_in ), file=sys.stderr )
        sys.exit( -1 )
    index1 = row1_index
    index2 = row2_index + 11
    for i in range( index1, index2 + 1 ):
      index_list.append( i )
  return( check_index_list( index_list, element_coordinates ) )


def parse_columns( string_in, element_coordinates = [ None, None ] ):
  """
  Convert a column specification to an index string.
  Acceptable column specifications include
    o  single column
         o  5
    o  column range
         o  6-8
    o  single and/or ranges of column separated by commas
         o  9-11,3
  """
  index_list = []
  string_in = re.sub( r'\s', '', string_in )
  for col_range in string_in.split( ',' ):
    mobj = re.match( r'([1-9][0-2]?)([-:]([1-9][0-2]?))?$', col_range )
    if( not mobj ):
      print( 'Error: spreadsheet cell: %s%s: bad column or column range \'%s\'' % ( element_coordinates[0], element_coordinates[1], col_range ), file=sys.stderr )
      sys.exit( -1 )

    col1 = int( mobj.group( 1 ) )
    col1_index = well_to_index( 1, 'A', col1, False, element_coordinates )
    col2_index = col1_index
    if( col1 < 1 or col1 > 12 ):
      print( 'Error: spreadsheet cell: %s%s: bad column value: \'%d\'' % ( element_coordinates[0], element_coordinates[1], col1 ), file=sys.stderr )
      sys.exit( -1 )
    if( mobj.group( 2 ) ):
      col2 = int( mobj.group( 3 ) )
      if( col2 < 1 or col2 > 12 ):
        print( 'Error: spreadsheet cell: %s%s: bad column value: \'%d\'' % ( element_coordinates[0], element_coordinates[1], col1 ), file=sys.stderr )
        sys.exit( -1 )
      if( col2 < col1 ):
        print( 'Error: spreadsheet cell: %s%s: bad column range: \'%s\'' % ( element_coordinates[0], element_coordinates[1], string_in ), file=sys.stderr )
        sys.exit( -1 )
      col2_index = well_to_index( 1, 'A', col2, False, element_coordinates )
    index1 = col1_index
    index2 = col2_index + 7
    for i in range( index1, index2 + 1 ):
      index_list.append( i )
  return( check_index_list( index_list, element_coordinates ) )


def check_rows( column_name_list, csv_rows ):
  """
  Trim off empty cells at end of rows and columns
  and check for internal empty cells. Allow empty
  internal row.
  Notes:
    o  we expect
         o  nrows 
  """
  # check for internal empty row
  csv_rows_out = []
  num_col = len( column_name_list )
  for irow, row_elements in enumerate( csv_rows ):
    if(len(row_elements) < num_col):
      print('Error: missing %d cell(s) in row %d' % ( num_col - len(row_elements), irow))
      sys.exit( -1 )
    # Number of cells that are supposed to have
    # content but don't.
    num_empty = 0
    row_elements_out = []
    for icol, cell in enumerate( row_elements ):
      if( icol == num_col ):
        break
      # Add valid cell to output row.
      if( len( cell ) > 0 or ( column_name_list[icol]['type'] in column_allow_empty_cell ) ):
        row_elements_out.append( cell )
      else:
        num_empty += 1
    # Save the row if there are no invalid empty cells.
    if( num_empty == 0 ):
      csv_rows_out.append( row_elements_out )
    elif( num_empty > 0 and num_empty < num_col ):
      # Note: we allow for an empty row with
      #       num_empty == num_col.
      print( 'Error: row %d has empty cells' % ( irow + 1 ) )
      sys.exit( -1 )
  return( csv_rows_out )


def read_samplesheet( file ):
  """
  Read CSV samplesheet input file.
  Notes:
    o  the first row in the file must have column header names.
    o  the column header names must be in the list 'column_header_name_list'.
    o  the column order is arbitrary.
  """
  samplesheet_row_list = []
  csv_rows = csv.reader( file, delimiter=',', quotechar='"')
  csv_rows = list( csv_rows )
  row_header = csv_rows[0]
  column_name_list = parse_header( row_header )
  csv_rows = check_rows( column_name_list, csv_rows )
  for row_elements in csv_rows[1:]:
    samplesheet_row_list.append( row_elements )
  return( column_name_list, samplesheet_row_list )


def check_sample_names( column_name_list, samplesheet_row_list ):
  """
  Check for name degeneracy.
  Check sample names for unacceptable characters and, if present, convert them to '.'.
  Sample names must begin with [a-zA-Z].
  Check that peak groups are positive integers.
  Unacceptable characters are characters that are not a-z, A-Z, 0-9, and '.'
  Check for name degeneracy after substitutions.
  Check that the barnyard sample is labeled 'Barnyard'.
  """
  sample_name_in_dict = {}
  sample_name_out_dict = {}
  num_sample_name = 0
  for row_elements in samplesheet_row_list:
    for i in range( len( row_elements ) ):
      column_name_dict = column_name_list[i]
      element_string = row_elements[i]
      if( column_name_dict['type'] != 'sample_name' ):
        continue
      sample_name_in_dict.setdefault( element_string, 0 )
      sample_name_in_dict[element_string] += 1
      num_sample_name += 1
      mobj = re.match( r'[a-zA-Z0-9]', element_string )
      if( not mobj ):
        print( 'Error: sample names must begin with an alphanumeric character', file=sys.stderr )
        sys.exit( -1 )
      # Convert invalid characters to periods.
      row_elements[i] = re.sub( r'[^a-zA-Z0-9.]', '.', element_string )
      sample_name_out_dict.setdefault( element_string, True )
  errorFlag = False
  for sample_name in sample_name_in_dict.keys():
    if( sample_name_in_dict[sample_name] > 1 ):
      print( 'Warning: sample name \'%s\' not unique. It is used %d times.' % ( sample_name, sample_name_in_dict[sample_name] ), file=sys.stderr )
  if( len( sample_name_out_dict ) != len( sample_name_in_dict ) ):
    print( 'Error: unacceptable names are not distinct after editing', file=sys.stderr )
    errorFlag = True
  if( errorFlag ):
    sys.exit( -1 )
  # Check barnyard sample name. (This is unnecessary, I believe.)
#   for row_elements in samplesheet_row_list:
#     for i in range( len( row_elements ) ):
#       column_name_dict = column_name_list[i]
#       element_string = row_elements[i]
#       if( column_name_dict['type'] != 'sample_name' ):
#         continue
#       mobj = re.search( r'barn', element_string.lower() )
#       if( mobj and element_string != 'Barnyard' ):
#         print( '**' )
#         print( '** Warning: barnyard sample name (%s) not \'Barnyard\'.' % ( element_string ), file=sys.stderr )
#         print( '**          Consider re-naming it to \'Barnyard\' for compatibility' )
#         print( '**          with the experiment dashboard.' )
#         print( '**' )
#       break
  return( samplesheet_row_list )


def check_genome_names( column_name_list, samplesheet_row_list ):
  """
  Check genome names and warn if not in our list.
  """
  missing_genome_names_dict = {}
  errorFlag = False
  for irow, row_elements in enumerate(samplesheet_row_list):
    for i in range( len( row_elements ) ):
      column_name_dict = column_name_list[i]
      element_string = row_elements[i]
      if( column_name_dict['type'] != 'genome' ):
        continue
      if( has_space(element_string) ):
        print( 'Error: genome name \'%s\' in row %d has one or more spaces' % ( element_string, irow ) )
        errorFlag = 1
      if( not row_elements[i] in genome_name_list ):
        missing_genome_names_dict.setdefault( row_elements[i], True )
  if( len( missing_genome_names_dict.keys() ) > 0 ):
    print( 'The following genomes are not in my list of known genomes (they may be mis-spelled or not in my list).', file=sys.stderr )
    for missing_genome_name in missing_genome_names_dict.keys():
      print( '  \'%s\'' % ( missing_genome_name ), file=sys.stderr )
  if( errorFlag ):
    sys.exit( -1 )
  return( 0 )


def check_peak_groups( column_name_list, samplesheet_row_list ):
  """
  Check peak group names and exit on error.
  """
  bad_peak_groups_dict = {}
  for row_elements in samplesheet_row_list:
    for i in range( len( row_elements ) ):
      column_name_dict = column_name_list[i]
      element_string = row_elements[i]
      if( column_name_dict['type'] != 'peak_group' ):
        continue
      if( len( element_string ) > 0 and re.search(r'[^a-zA-Z0-9_]', element_string ) ):
        bad_peak_groups_dict.setdefault( element_string, True )
  if( len( bad_peak_groups_dict.keys() ) > 0 ):
    print('Unacceptable peak group names (must use only alphabetic, positive integer, and underscore characters):')
    for bad_peak_group in bad_peak_groups_dict.keys():
      print( '  \'%s\'' % ( bad_peak_group ) )
    sys.exit( -1 )
  return( 0 )


def check_peak_files( column_name_list, samplesheet_row_list ):
  """
  Check peak file names and exit on error.
  """
  bad_peak_files_dict = {}
  for row_elements in samplesheet_row_list:
    for i in range( len( row_elements ) ):
      column_name_dict = column_name_list[i]
      element_string = row_elements[i]
      if( column_name_dict['type'] != 'peak_file' ):
        continue
      if( ( len( element_string ) > 0 ) and ( re.match(r'[^/]', element_string ) or has_space( element_string ) ) ):
        bad_peak_files_dict.setdefault( element_string, True )
  if( len( bad_peak_files_dict.keys() ) > 0 ):
    print('Unacceptable peak file names (must start with \'/\' and must not contain spaces):')
    for bad_peak_file in bad_peak_files_dict.keys():
      print( '  \'%s\'' % ( bad_peak_file ) )
    sys.exit( -1 )
  return( 0 )


def check_peak_spec( column_name_list, samplesheet_row_list ):
  """
  Check that each sample has a peak_group or a peak_file or both.
  """
  bad_peak_dict = {}
  for row_elements in samplesheet_row_list:
    peak_group_flag = False
    peak_file_flag = False
    for i in range( len( row_elements ) ):
      column_name_dict = column_name_list[i]
      element_string = row_elements[i]
      if( column_name_dict['type'] == 'sample_name' ):
        sample_name = element_string
      if( column_name_dict['type'] == 'peak_group' and len( element_string ) > 0 ):
        peak_group_flag = True
      if( column_name_dict['type'] == 'peak_file' and len( element_string ) > 0 ):
        peak_file_flag = True
    if( ( not peak_group_flag ) and ( not peak_file_flag ) ):
      bad_peak_dict.setdefault( sample_name, True )
  if( len( bad_peak_dict.keys() ) > 0 ):
    print('Samples have no peak_group and no peak_file values:')
    for bad_peak_spec in bad_peak_dict.keys():
      print( '  \'%s\'' % ( bad_peak_spec ) )
    sys.exit( -1 )
  return( 0 )


def check_external_sample_name( column_name_list, samplesheet_row_list ):
  """
  Check that each external_sample_name is valid.
  """
  errorFlag = False
  sample_name_in_dict = {}
  sample_name_out_dict = {}
  for irow, row_elements in enumerate(samplesheet_row_list):
    for i in range( len( row_elements ) ):
      column_name_dict = column_name_list[i]
      element_string = row_elements[i]
      if( column_name_dict['type'] != 'external_sample_name' ):
        continue
      # Allow empty external sample name cell.
      if( len( element_string ) == 0 ):
        continue
      sample_name_in_dict.setdefault( element_string, 0 )
      sample_name_in_dict[element_string] += 1
      mobj = re.match( r'[a-zA-Z0-9]', element_string )
      if( not mobj ):
        print( 'Error: bad external sample name in row %d (must begin with an alphanumeric character)' % ( irow ), file=sys.stderr )
        errorFlag = True
      sample_name_out_dict.setdefault( element_string, True )
  for sample_name in sample_name_in_dict.keys():
    if( sample_name_in_dict[sample_name] > 1 ):
      print( 'Warning: external sample name \'%s\' not unique. It is used %d times.' % ( sample_name, sample_name_in_dict[sample_name] ), file=sys.stderr )
  if( len( sample_name_out_dict ) != len( sample_name_in_dict ) ):
    print( 'Error: unacceptable names are not distinct after editing', file=sys.stderr )
    errorFlag = True
  if( errorFlag ):
    sys.exit( -1 )
  return( samplesheet_row_list )


def check_tissue( column_name_list, samplesheet_row_list ):
  """
  Check that sample tissue strings are valid.
  There are no additional restrictions on tissue
  names at this time so this function has no
  function.
  """
  bad_tissue_dict = {}
  for row_elements in samplesheet_row_list:
    for i in range( len( row_elements ) ):
      column_name_dict = column_name_list[i]
      element_string = row_elements[i]
      if( column_name_dict['type'] != 'tissue' ):
        continue
      # Allow empty external sample name cell.
      if( len( element_string ) == 0 ):
        continue
  if( len( bad_tissue_dict.keys() ) > 0 ):
    print('Error: invalid tissue names:')
    for bad_tissue in bad_tissue_dict.keys():
      print( '  \'%s\'' % ( bad_tissue ) )
    sys.exit( -1 )
  return( samplesheet_row_list )


def check_wrap_group( column_name_list, samplesheet_row_list ):
  """
  Check that wrap_group values are valid.
  """
  bad_wrap_group_dict = {}
  for row_elements in samplesheet_row_list:
    wrap_group_flag = False
    for i in range( len( row_elements ) ):
      column_name_dict = column_name_list[i]
      element_string = row_elements[i]
      if( column_name_dict['type'] != 'wrap_group' ):
        continue
      sample_name = element_string
      if( len( element_string ) > 0 and re.search(r'[^-_.a-zA-Z0-9]', element_string ) ):
        bad_wrap_group_dict.setdefault( element_string, True )
  if( len( bad_wrap_group_dict.keys() ) > 0 ):
    print('Sample wrap_group values have invalid characters.' )
    print('Valid characters are alphabetic, positive integers,')
    print('and ".", "_", and "-".')
    for bad_wrap_group in bad_wrap_group_dict.keys():
      print( '  \'%s\'' % ( bad_wrap_group ) )
    sys.exit( -1 )
  return( 0 )


def get_peak_files(  column_name_list, samplesheet_row_list ):
  peak_files_dict = {}

  for row_elements in samplesheet_row_list:
    sample_name = None
    peak_file = None
    for i in range( len( row_elements ) ):
      column_name_dict = column_name_list[i]
      element_string = row_elements[i]
      if( column_name_dict['type'] == 'sample_name' ):
        sample_name = element_string
      if( column_name_dict['type'] == 'peak_file' ):
        peak_file = element_string
    if( sample_name != None and peak_file != None ):
      if( peak_files_dict.get( 'sample_name' ) != None and
          peak_file != peak_files_dict['peak_file'] ):
        print( 'Error: inconsistent peak file names for sample: %s' % (sample_name) )
        sys.exit( -1 )
      else:
        peak_files_dict[sample_name] = peak_file
  return( peak_files_dict )


def expand_rows( string_in, element_coordinates = [ None, None ] ):
  """
  Expand a P7 row specification to a list of rows.
  Acceptable row specifications include
    o  single row
         o  B
    o  row range
         o  E-G
    o  single and/or ranges of rows separated by commas
         o  E-F,H
  """
  string_in = re.sub( r'\s', '', string_in )
  row_list = []
  for row_range in string_in.split( ',' ):
    mobj = re.match( r'([a-hA-H])([-:]([a-hA-H]))?$', row_range )
    if( not mobj ):
      print( 'Error: spreadsheet cell: %s%s: bad row or row range \'%s\'' % ( element_coordinates[0], element_coordinates[1], row_range ), file=sys.stderr )
      sys.exit( -1 )
    icode_row1 = ord(mobj.group( 1 ))
    icode_row2 = icode_row1
    if( mobj.group( 2 ) ):
      icode_row2 = ord(mobj.group( 3 ))
      if( icode_row2 < icode_row1 ):
        print( 'Error: spreadsheet cell: %s%s: bad row range: \'%\'' % ( element_coordinates[0], element_coordinates[1], string_in ), file=sys.stderr )
        sys.exit( -1 )
    for icode in range(icode_row1, icode_row2+1):
      row_list.append(chr(icode))
  return(row_list)


def expand_columns( string_in, element_coordinates = [ None, None ] ):
  """
  Expand a P5 column specification to a list of columns.
  Acceptable column specifications include
    o  single column
         o  5
    o  column range
         o  6-8
    o  single and/or ranges of column separated by commas
         o  9-11,3
  """
  string_in = re.sub( r'\s', '', string_in )
  column_list = []
  for col_range in string_in.split( ',' ):
    mobj = re.match( r'([1-9][0-2]?)([-:]([1-9][0-2]?))?$', col_range )
    if( not mobj ):
      print( 'Error: spreadsheet cell: %s%s: bad column or column range \'%s\'' % ( element_coordinates[0], element_coordinates[1], col_range ), file=sys.stderr )
      sys.exit( -1 )

    icol_col1 = int( mobj.group( 1 ) )
    icol_col2 = icol_col1
    if( icol_col1 < 1 or icol_col1 > 12 ):
      print( 'Error: spreadsheet cell: %s%s: bad column value: \'%d\'' % ( element_coordinates[0], element_coordinates[1], col1 ), file=sys.stderr )
      sys.exit( -1 )
    if( mobj.group( 2 ) ):
      icol_col2 = int( mobj.group( 3 ) )
      if( icol_col2 < 1 or icol_col2 > 12 ):
        print( 'Error: spreadsheet cell: %s%s: bad column value: \'%d\'' % ( element_coordinates[0], element_coordinates[1], col1 ), file=sys.stderr )
        sys.exit( -1 )
      if( icol_col2 < icol_col1 ):
        print( 'Error: spreadsheet cell: %s%s: bad column range: \'%s\'' % ( element_coordinates[0], element_coordinates[1], string_in ), file=sys.stderr )
        sys.exit( -1 )
    for i in range( icol_col1, icol_col2 + 1 ):
      column_list.append(str(i))
  return(column_list)


def test_pcr_format(column_name_list):
  """
  Find the P7 and P5 columns in the header column_name_list.
  Return the column indices, or None if the P7 and P5 are
  not in the row and column format.
  """
  pcr7_column = None
  pcr5_column = None
  for icol in range(len(column_name_list)):
    column_dict = column_name_list[icol]
    if(column_dict['type'] == 'p7' and
       column_dict['format'] == 'rows'):
      pcr7_column = icol
    elif(column_dict['type'] == 'p5' and
         column_dict['format'] == 'columns'):
      pcr5_column = icol
  return(pcr7_column, pcr5_column)


def expand_sample_rows(column_name_list, samplesheet_row_list):
  """
  Expand sample rows by each P7 row and P5 column listed for
  each sample. For example, if sample_A has P7 rows 'A,D'
  and P5 columns '5,3', expand_sample_rows returns a list
  that includes two rows for sample_A where the first sample
  row has A for the P7 row and 5 for the P5 column and the
  second row has D and 3 for the P7 row and P5 column,
  respectively.
  """
  # If P7 format != p7_rows or P5 format != p5_columns, then
  # return without expanding.
  pcr7_column, pcr5_column = test_pcr_format(column_name_list)
  if(pcr7_column == None or pcr5_column == None):
    return(samplesheet_row_list)

  # Expand samples by pairing the i-th P7 row with the i-th P5
  # column for each of the specified rows and columns. Each sample
  # row has one PCR row and one PCR column.
  new_samplesheet_row_list = []
  num_element = len(column_name_list)
  for irow, samplesheet_row in enumerate(samplesheet_row_list):
    element_coordinates = [str( irow + 2 ), chr( pcr7_column + ord( 'A' ))]
    row_list = expand_rows(samplesheet_row[pcr7_column], element_coordinates = element_coordinates)
    element_coordinates = [str( irow + 2 ), chr( pcr5_column + ord( 'A' ))]
    column_list = expand_columns(samplesheet_row[pcr5_column], element_coordinates = element_coordinates)
    if(len(row_list) != len(column_list)):
      print('Error: number of rows is not equal to number of columns in sample sheet row %d' % (irow + 2))
      sys.exit(-1)
    num_expand = len(row_list)
    for iexpand in range(num_expand):
      element_list = []
      for ielem in range(num_element):
        if(ielem == pcr7_column):
          element_list.append(row_list[iexpand])
        elif(ielem == pcr5_column):
          element_list.append(column_list[iexpand])
        else:
          element_list.append(samplesheet_row[ielem])
      new_samplesheet_row_list.append(element_list)
  return(new_samplesheet_row_list)


def make_samplesheet_indexes( column_name_list, samplesheet_row_list ):
  """
  Make well index lists for N5, N7, P5, P7 barcode wells from the input samplesheet information.
  """
  num_col = len( column_name_list )
  row_out_list = []
  for irow, row_elements in enumerate( samplesheet_row_list ):
    if( len( row_elements ) < num_col ):
      print( 'Error: missing cells in row %d: %s' % ( irow + 1, ', '.join('"{0}"'.format(e) for e in row_elements ) ), file=sys.stderr )
      sys.exit( -1 )
    icol = 0
    # Initialize optional column values.
    wrap_group = None
    tissue = None
    for element_string, column_name_dict in zip( row_elements, column_name_list ):
      icol += 1

      element_coordinates = [ str( irow + 2 ), chr( icol + ord( 'A' ) - 1 ) ]
      if( column_name_dict['type'] == 'n7' ):
        max_index = 384
        if( column_name_dict['format'] == 'indexes' ):
          n7_index_list = parse_indexes( element_string, max_index, element_coordinates )
        elif( column_name_dict['format'] == 'wells' ):
          n7_index_list = parse_wells( element_string, True, max_index, element_coordinates )
        else:
          print( 'Error: unexpected N7 format', file=sys.stderr )
          sys.exit( -1 )
      elif( column_name_dict['type'] == 'n5' ):
        max_index = 384
        if( column_name_dict['format'] == 'indexes' ):
          n5_index_list = parse_indexes( element_string, max_index, element_coordinates )
        elif( column_name_dict['format'] == 'wells' ):
          n5_index_list = parse_wells( element_string, False, max_index, element_coordinates )
        else:
          print( 'Error: unexpected N5 format', file=sys.stderr )
          sys.exit( -1 )
      elif( column_name_dict['type'] == 'p7' ):
        max_index = 96
        if( column_name_dict['format'] == 'indexes' ):
          p7_index_list = parse_indexes( element_string, max_index, element_coordinates )
        elif( column_name_dict['format'] == 'wells' ):
          p7_index_list = parse_wells( element_string, True, max_index, element_coordinates )
        elif( column_name_dict['format'] == 'rows' ):
          p7_index_list = parse_rows( element_string, element_coordinates )
        else:
          print( 'Error: unexpected P7 format', file=sys.stderr )
          sys.exit( -1 )
      elif( column_name_dict['type'] == 'p5' ):
        max_index = 96
        if( column_name_dict['format'] == 'indexes' ):
          p5_index_list = parse_indexes( element_string, max_index, element_coordinates )
        elif( column_name_dict['format'] == 'wells' ):
          p5_index_list = parse_wells( element_string, False, max_index, element_coordinates )
        elif( column_name_dict['format'] == 'columns' ):
          p5_index_list = parse_columns( element_string, element_coordinates )
        else:
          print( 'Error: unexpected P5 format', file=sys.stderr )
          sys.exit( -1 )
      elif( column_name_dict['type'] == 'sample_name' ):
          sample_name = element_string
      elif( column_name_dict['type'] == 'genome' ):
          genome = element_string
      elif( column_name_dict['type'] == 'peak_group' ):
          peak_group = element_string
      elif( column_name_dict['type'] == 'peak_file' ):
          peak_file = element_string
      elif( column_name_dict['type'] == 'external_sample_name' ):
          external_sample_name = element_string
      elif( column_name_dict['type'] == 'tissue' ):
          tissue = element_string
      elif( column_name_dict['type'] == 'wrap_group' ):
          wrap_group = element_string
    #
    row_out_list.append( { 'sample_name': sample_name,
                           'n7_index_list': n7_index_list,
                           'p7_index_list': p7_index_list,
                           'n5_index_list': n5_index_list,
                           'p5_index_list': p5_index_list,
                           'genome': genome,
                           'peak_group': peak_group,
                           'peak_file': peak_file,
                           'external_sample_name': external_sample_name,
                           'tissue': tissue,
                           'wrap_group': wrap_group } )
  return( row_out_list )


def check_sample_identifier( row_out_list, sample_identifier ):
  """
  Check for samples that have in common sample identifier well indexes.
  Typically, sci-ATAC experiments use the N5 barcode to identify samples. Check whether
  these sample identifier barcodes occur in more than one sample.
  Notes:
    o  some experiments use a combination of ligation and PCR barcodes to
       identify samples, and the same ligation barcodes can be assigned
       intentionally to more than one sample 
    o  issue warning and continue rather than error and exit
  """
  key = '%s_index_list' % ( sample_identifier )
  index_dict = {}
  for row_out in row_out_list:
    for i in row_out[key]:
      index_dict.setdefault( i, 0 )
      index_dict[i] += 1
  duplicate_list = []
  for i in index_dict.keys():
    if( index_dict[i] > 1 ):
      duplicate_list.append( i )

  if( len( duplicate_list ) > 0 ):
    print( 'Warning: the following samples have %s wells in common' % ( sample_identifier ), file=sys.stderr )
    for row_out in row_out_list:
      if( set( duplicate_list ).intersection( set( row_out[key] ) ) ):
        print( '  %s' % ( row_out['sample_name'] ), file=sys.stderr )
  return( 0 )


def dump_row_out_list( row_out_list ):
  """
  Diagnostic function to dump (barcode) well index lists.
  """
  for row_out in row_out_list:
    print( 'sample_name: %s' % ( row_out['sample_name'] ) )
    print( '  n7_index_list: %s' % ( make_index_string( row_out['n7_index_list'] ) ) )
    print( '  p7_index_list: %s' % ( make_index_string( row_out['p7_index_list'] ) ) )
    print( '  p5_index_list: %s' % ( make_index_string( row_out['p5_index_list'] ) ) )
    print( '  n5_index_list: %s' % ( make_index_string( row_out['n5_index_list'] ) ) )
    print( '  genome: %s' % ( row_out['genome'] ) )
  return( 0 )


def write_samplesheet_index_format( file, row_out_list ):
  """
  Write a well index samplesheet file in Andrew's format.
  """
  print( 'sample_id\tranges\tgenome', file=file )
  for row_out in row_out_list:
    print( '%s\t%s:%s:%s:%s\t%s' % ( row_out['sample_name'],
                                     make_index_string( row_out['n7_index_list'] ),
                                     make_index_string( row_out['p7_index_list'] ),
                                     make_index_string( row_out['p5_index_list'] ),
                                     make_index_string( row_out['n5_index_list'] ),
                                     row_out['genome']), file=file )

  return( 0 )


def get_pcr_row_col( column_name_list, samplesheet_row_list ):
  """
  Given lists of PCR rows and columns by sample, return lists of
  distinct PCR rows and columns for JSON output file. This values
  are used by the demux_dash.
  """
  p5_re_pattern = r'([1-9][0-2]?)([-:]([1-9][0-2]?))?$'
  p7_re_pattern = r'([a-hA-H])([-:]([a-hA-H]))?$'

  # find required samplesheet column
  p5_samplesheet_col = None
  p7_samplesheet_col = None
  for icol, column_name_dict in enumerate( column_name_list ):
    if( column_name_dict['type'] == 'p5' ):
      p5_samplesheet_col = icol
    if( column_name_dict['type'] == 'p7' ):
      p7_samplesheet_col = icol

  # gather plate values
  pair_samplesheet_row_list = []
  for samplesheet_row in samplesheet_row_list:
    p5_samplesheet_row_list = []
    for value_range in samplesheet_row[p5_samplesheet_col].split( ',' ):
      value_range = re.sub( r'\s', '', value_range )
      mobj = re.match( p5_re_pattern, value_range )
      if( not mobj ):
        print( 'Error: bad value or value range \'%s\'' % ( value_range ), file=sys.stderr )
        sys.exit( -1 )
      col1 = int( mobj.group( 1 ) )
      col2 = col1
      if( mobj.group( 2 ) ):
        col2 = int( mobj.group( 3 ) )
      for col in range( col1, col2 + 1 ):
        p5_samplesheet_row_list.append( str( col ) )

    p7_samplesheet_row_list = []
    for value_range in samplesheet_row[p7_samplesheet_col].split( ',' ):
      value_range = re.sub( r'\s', '', value_range )
      mobj = re.match( p7_re_pattern, value_range )
      if( not mobj ):
        print( 'Error: bad value or value range \'%s\'' % ( value_range ), file=sys.stderr )
        sys.exit( -1 )
      row1 = mobj.group( 1 )
      row2 = row1
      if( mobj.group( 2 ) ):
        row2 = mobj.group( 3 )
      for orow in range( ord( row1 ), ord( row2 ) + 1 ):
        p7_samplesheet_row_list.append( chr( orow ) )

    pair_samplesheet_row_list.append( [ p5_samplesheet_row_list, p7_samplesheet_row_list ] )

  pair_list_dict = {}
  for pair_samplesheet_row in pair_samplesheet_row_list:
    key = '_'.join( pair_samplesheet_row[0] ) + '_' + '_'.join( pair_samplesheet_row[1] )
    pair_list_dict.setdefault( key, pair_samplesheet_row )

  p5_col_list = []
  p7_row_list = []
  for key in pair_list_dict.keys():
    p5_col_list.extend( pair_list_dict[key][0] )
    p7_row_list.extend( pair_list_dict[key][1] )

  return( p5_col_list, p7_row_list )


def get_pcr_wells( column_name_list, samplesheet_row_list ):
  """
  Given the samplesheet rows, return lists of specified PCR wells,
  if any.
  """
  p5_samplesheet_col = None
  p7_samplesheet_col = None
  for icol, column_name_dict in enumerate( column_name_list ):
    if( column_name_dict['type'] == 'p5' ):
      p5_samplesheet_col = icol
    if( column_name_dict['type'] == 'p7' ):
      p7_samplesheet_col = icol

  pair_samplesheet_row_list = []
  for samplesheet_row in samplesheet_row_list:
    p5_samplesheet_row_list = []
    for value_range in samplesheet_row[p5_samplesheet_col].split( ',' ):
      value_range = re.sub( r'\s', '', value_range )
      mobj = re.match( r'([pP][0]?([1-4*])[-:])?([a-hA-H])([0]?[1-9][0-2]?)([-:]([pP][0]?([1-4*])[-:])?([a-hA-H])([0]?[1-9][0-2]?))?$', value_range )
      if( not mobj ):
        print( 'Error: bad well or well range \'%s\'' % ( value_range ), file=sys.stderr )
        sys.exit( -1 )

      # Assume that the well-range strings were checked earlier in parse_wells() call so
      # those tests are omitted here. Also, this is PCR for which there is a single plate=1.
      # first well
      row1 = mobj.group( 3 )
      col1 = int( mobj.group( 4 ) )
      row2 = row1
      col2 = col1
      if( col1 < 1 or col1 > 12 ):
        print( 'Error: bad well: \'%s\'' % ( value_range ), file=sys.stderr )
        sys.exit( -1 )
      if( mobj.group( 5 ) ):
        # second well, if this is a range
        row2 = mobj.group( 8 )
        col2 = int( mobj.group( 9 ) )
      across_row_first = False # p5
      index1 = well_to_index( 1, row1, col1, across_row_first, [ None, None ] )
      index2 = well_to_index( 1, row2, col2, across_row_first, [ None, None ] )

      for well_index in range( index1, index2 + 1 ):
        well = index_to_well( well_index - 1, across_row_first )
        p5_samplesheet_row_list.append( well )

    p7_samplesheet_row_list = []
    for value_range in samplesheet_row[p7_samplesheet_col].split( ',' ):
      value_range = re.sub( r'\s', '', value_range )
      mobj = re.match( r'([pP][0]?([1-4*])[-:])?([a-hA-H])([0]?[1-9][0-2]?)([-:]([pP][0]?([1-4*])[-:])?([a-hA-H])([0]?[1-9][0-2]?))?$', value_range )
      if( not mobj ):
        print( 'Error: bad well or well range \'%s\'' % ( value_range ), file=sys.stderr )
        sys.exit( -1 )

      # Assume that the well-range strings were checked earlier in parse_wells() call so
      # those tests are omitted here. Also, this is PCR for which there is a single plate=1.
      # first well
      row1 = mobj.group( 3 )
      col1 = int( mobj.group( 4 ) )
      row2 = row1
      col2 = col1
      if( col1 < 1 or col1 > 12 ):
        print( 'Error: bad well: \'%s\'' % ( value_range ), file=sys.stderr )
        sys.exit( -1 )
      if( mobj.group( 5 ) ):
        # second well, if this is a range
        row2 = mobj.group( 8 )
        col2 = int( mobj.group( 9 ) )
      across_row_first = True # p7
      index1 = well_to_index( 1, row1, col1, across_row_first, [ None, None ] )
      index2 = well_to_index( 1, row2, col2, across_row_first, [ None, None ] )

      for well_index in range( index1, index2 + 1 ):
        well = index_to_well( well_index - 1, across_row_first )
        p7_samplesheet_row_list.append( well )

    pair_samplesheet_row_list.append( [ p5_samplesheet_row_list, p7_samplesheet_row_list ] )

  pair_list_dict = {}
  for pair_samplesheet_row in pair_samplesheet_row_list:
    key = '_'.join( pair_samplesheet_row[0] ) + '_' + '_'.join( pair_samplesheet_row[1] )
    pair_list_dict.setdefault( key, pair_samplesheet_row )

  p5_well_list = []
  p7_well_list = []
  for key in pair_list_dict.keys():
    p5_well_list.extend( pair_list_dict[key][0] )
    p7_well_list.extend( pair_list_dict[key][1] )

  return( p5_well_list, p7_well_list )


def write_samplesheet_json_format( file, column_name_list, samplesheet_row_list, row_out_list, wrap_groups_dict, level = 3, number_wells = 384, tn5_barcodes = False, use_all_barcodes = False, illumina_run_directory = 'NA', hash_file = None ):
  """
  Write an output samplesheet file in JSON format.
  """
  # Store input samplesheet for for reference if questions arise.
  input_samplesheet_rows = []
  for samplesheet_row in samplesheet_row_list:
    input_samplesheet_rows.append( ','.join( '"{0}"'.format( e ) for e in samplesheet_row ) )

  # Store sample information for processing pipeline.
  sample_index_list = []
  for row_out in row_out_list:
    sample_index_list.append( { 'sample_id' : row_out['sample_name'],
                                'ranges' : ':'.join( [ make_index_string( row_out['n7_index_list'] ),
                                                       make_index_string( row_out['p7_index_list'] ),
                                                       make_index_string( row_out['p5_index_list'] ),
                                                       make_index_string( row_out['n5_index_list'] ) ] ),
                                'genome' : row_out['genome'],
                                'peak_group' : row_out['peak_group'],
                                'peak_file' : row_out['peak_file'] })

  # Store information for dashboard(s).

  # Note: assume that the header was checked for consistent
  #       p5 => columns and p7 => rows
  pcr_format = None
  for icol, column_name_dict in enumerate( column_name_list ):
    if( column_name_dict['type'] == 'p5' ):
      if( column_name_dict['format'] == 'columns' ):
        pcr_format = 'row_col'
      elif( column_name_dict['format'] == 'indexes' ):
        pcr_format = 'indexes'
      elif( column_name_dict['format'] == 'wells' ):
        pcr_format = 'wells'

  # PCR rows and columns specified?
  p5_col_list = None
  p7_row_list = None
  if( pcr_format == 'row_col' ):
    p5_col_list, p7_row_list = get_pcr_row_col( column_name_list, samplesheet_row_list )

  # PCR wells specified?
  p5_well_list = None
  p7_well_list = None
  if( pcr_format == 'wells' ):
    p5_well_list, p7_well_list = get_pcr_wells( column_name_list, samplesheet_row_list ) 

  # Tissue dict.
  tissue_dict = {}
  for row_out in row_out_list:
    tissue_dict[row_out['sample_name']] = row_out.get('tissue', '')
    
  # external_sample_name dict.
  external_sample_name_dict = {}
  for row_out in row_out_list:
    external_sample_name_dict[row_out['sample_name']] = row_out.get('external_sample_name', '')

#  # Make input column name list.
#  input_column_name_list = []
#  for column_name_dict in column_name_list:
#    if(column_name_dict.get('format', None) == None):
#      input_column_name_list.append(column_name_dict['type'])
#    else:
#      input_column_name_list.append(column_name_dict['type'] + '_' + column_name_dict['format'])
#  input_column_names = ','.join(input_column_name_list)

  # Make input column name list.
  input_column_name_list = []
  for column_name_dict in column_name_list:
    input_column_name_list.append(column_name_dict)
 


  # JSON structure.
  sample_data = { 'json_file_version' : json_file_version,
                  'illumina_run_directory' : illumina_run_directory,
                  'level' : level,
                  'number_wells' : number_wells,
                  'tn5_barcodes' : tn5_barcodes,
                  'use_all_barcodes' : use_all_barcodes,
                  'hash_file' : hash_file,
                  'input_samplesheet_column_names' : input_column_name_list,
                  'input_samplesheet_rows' : input_samplesheet_rows,
                  'pcr_format': pcr_format,
                  'p5_col_list' : p5_col_list,
                  'p7_row_list' : p7_row_list,
                  'p5_well_list' : p5_well_list,
                  'p7_well_list' : p7_well_list,
                  'sample_index_list' : sample_index_list,
                  'external_sample_name_dict': external_sample_name_dict,
                  'tissue_dict': tissue_dict,
                  'wrap_groups_dict': wrap_groups_dict,
                }
# add tissue list
# sample_name to external_sample_name dict

  file.write(json.dumps(sample_data, indent=4))

  return( 0 )


#
# Count distinct well indexes.
#
def count_wells( index_list ):
  num_well = len( set( index_list ) )
  return( num_well )


#
# Gather wrap group information. That is, the
# wrap group names, if they are given in the
# input samplesheet file, and the samples that belong
# to them.
#
def get_wrap_groups( row_out_list ):
  wrap_groups_dict = {}
  for irow, row_out in enumerate(row_out_list):
    # Gather wrap_group values in wrap_groups_dict.
    # Does this sample have a wrap group value?
    if( row_out.get( 'wrap_group' ) != None and row_out['wrap_group'] != '' ):
      if( wrap_groups_dict.get( row_out['wrap_group'] ) == None ):
        wrap_groups_dict[row_out['wrap_group']] = [ row_out['sample_name'] ]
      elif( row_out['sample_name'] not in wrap_groups_dict[row_out['wrap_group']] ):
        wrap_groups_dict[row_out['wrap_group']].append( row_out['sample_name'] )
  return( wrap_groups_dict )



def samplesheet_report( samplesheet_row_list, row_out_list, wrap_groups_dict, args ):
  print()
  print( '== Samplesheet information ==' )
  print( '  Tn5 barcodes:      %r' % ( args.tn5_barcodes ) )
  print( '  Level:             %s' % ( args.level ) )
  print( '  Number of wells:   %d' % ( args.number_wells ) )
  print( '  Sample identifier: %s' % ( args.sample_identifier ) )
  print( '  Use all barcodes:  %r' % ( args.use_all_barcodes ) )

  string_list = []
  print( '  Sample names after converting unacceptable characters to \'.\':' )
  for irow, row_out in enumerate(row_out_list):
    if( row_out_list[irow]['sample_name'] in string_list):
      continue
    print( '    %s' % ( row_out['sample_name'] ) )
    string_list.append(row_out_list[irow]['sample_name'])

  print( '  Sample well counts:' )
  max_len_samplename = 0
  for row_out in row_out_list:
    if( len( row_out['sample_name'] ) > max_len_samplename ):
      max_len_samplename = len( row_out['sample_name'] )
  if( len( 'name' ) > max_len_samplename ):
    max_len_samplename = len( 'name' )
  print( '    name%s    N7    P7    P5    N5' % ( ' ' * ( max_len_samplename - len( 'name' ) ) ) )
  for irow, row_out in enumerate(row_out_list):
    if(irow > 0
       and row_out_list[irow]['sample_name'] == row_out_list[irow-1]['sample_name']
       and count_wells(row_out_list[irow]['n7_index_list']) == count_wells(row_out_list[irow-1]['n7_index_list'])
       and count_wells(row_out_list[irow]['p7_index_list']) == count_wells(row_out_list[irow-1]['p7_index_list'])
       and count_wells(row_out_list[irow]['p5_index_list']) == count_wells(row_out_list[irow-1]['p5_index_list'])
       and count_wells(row_out_list[irow]['n5_index_list']) == count_wells(row_out_list[irow-1]['n5_index_list'])):
      continue

    print( '    %s%s    %d    %d    %d    %d' % ( row_out['sample_name'],
                                          ' ' * ( max_len_samplename - len( row_out['sample_name'] ) ),
                                          count_wells( row_out['n7_index_list'] ),
                                          count_wells( row_out['p7_index_list'] ),
                                          count_wells( row_out['p5_index_list'] ),
                                          count_wells( row_out['n5_index_list'] ) ) )
  print( '  Sample peak groups and files:' )
  max_len_peak_group = 0
  max_len_genome = 0
  for row_out in row_out_list:
    if( len( row_out['peak_group'] ) > max_len_peak_group ):
      max_len_peak_group = len( row_out['peak_group'] )
    if( len( row_out['genome'] ) > max_len_genome ):
      max_len_genome = len( row_out['genome'] )
  if( len( 'peak_group' ) > max_len_peak_group ):
    max_len_peak_group = len( 'peak_group' )
  if( len( 'genome' ) > max_len_genome ):
    max_len_genome = len( 'genome' )
  print( '    name%s    genome%s    peak_group%s    peak_file' % ( ' ' * ( max_len_samplename - len( 'name' ) ), ' ' * ( max_len_genome - len( 'genome' ) ), ' ' * ( max_len_peak_group - len( 'peak_group' ) ) ) )
  for irow, row_out in enumerate(row_out_list):
    if(irow > 0
       and row_out_list[irow]['sample_name'] == row_out_list[irow-1]['sample_name'] 
       and row_out_list[irow]['genome'] == row_out_list[irow-1]['genome']
       and row_out_list[irow]['peak_group'] == row_out_list[irow-1]['peak_group']
       and row_out_list[irow]['peak_file'] == row_out_list[irow-1]['peak_file']):
      continue

    print( '    %s%s    %s%s    %s%s    %s' % ( row_out['sample_name'],
                                          ' ' * ( max_len_samplename - len( row_out['sample_name'] ) ),
                                          row_out['genome'],
                                          ' ' * ( max_len_genome - len( row_out['genome'] ) ),
                                          row_out['peak_group'],
                                          ' ' * ( max_len_peak_group - len( row_out['peak_group'] ) ),
                                          row_out['peak_file'] ) )

  # Report external sample names.
  print( '  External sample names with tabs converted to spaces:' )
  for irow, row_out in enumerate(row_out_list):
    if(irow == 0 or row_out['sample_name'] != row_out_list[irow-1]['sample_name']):
      print('    %s\t%s' % ( row_out['sample_name'], row_out.get('external_sample_name', '')))

  # Report tissue names.
  print( '  Tissue assignments with tabs converted to spaces:' )
  for irow, row_out in enumerate(row_out_list):
    if(irow == 0 or row_out['sample_name'] != row_out_list[irow-1]['sample_name']):
      print('    %s\t%s' % ( row_out['sample_name'], row_out.get('tissue', '')))

  # Report information about wrap groups for the wrapping,
  # if it exists in the input samplesheet file.
  if( len( wrap_groups_dict.keys() ) > 0 ):
    print( '  Distribution groups:' )
    for wrap_group in wrap_groups_dict.keys():
      print( '   Group: %s' % ( wrap_group ) )
      print( '     Samples:' )
      for sample_name in wrap_groups_dict[wrap_group]:
        print( '        %s' % ( sample_name ) )

  print( '  Illumina run directory: %s' % ( args.run_dir ) )
  print( '  Run sciatac_samplesheet.py -d for more information.' )
  return( 0 )


def write_samplesheet_template():
  filename = 'samplesheet.template.csv'
  with open( filename, 'wt' ) as fp:
    print( 'n7_wells,p7_rows,p5_columns,n5_wells,sample_name,genome,peak_group,peak_file,external_sample_name,tissue,wrap_group', file=fp )
    print( 'p1:A01-p1:H12,a-d,5-8,p1:a01-p1:h01,sample1,barnyard,group_1,,,,None', file=fp )
  return( 0 )


if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='A program to convert sci-ATAC CSV samplesheet to pipeline samplesheet.')
  parser.add_argument('-i', '--input', required=False, default=None, help='Input CSV samplesheet filename (required string).')
  parser.add_argument('-o', '--output', required=False, default=None, help='Output samplesheet filename (required string).')
  parser.add_argument('-f', '--format', required=False, choices=[ 'json', 'index'], default='json', help='Output file format (default: \'%(default)s\') (optional string).')
  parser.add_argument('-r', '--run_dir', required=False, default=None, help='Illumina run directory path (optional string).')
  parser.add_argument('-l', '--level', type=int, required=False, choices=[ 2, 3 ], default=3, help='Two or three level sci-ATAC-seq experiment (default: %(default)s) (optional integer).')
  parser.add_argument('-w', '--number_wells', type=int, required=False, choices=[ 96, 384 ], default=384, help='Number of barcode set wells (default: %(default)s) (optional integer).')
  parser.add_argument('-t', '--tn5_barcodes', required=False, action='store_true', help='Tn5 has barcodes (optional flag).')
  parser.add_argument('-s', '--sample_identifier', required=False, choices=[ 'n5', 'n7' ], default='n5', help='Ligation barcode that identifies the sample (default: \'%(default)s\') used to check for duplicates (optional string).')
  parser.add_argument('--use_all_barcodes', required=False, action='store_true', help='Use all barcodes to demultiplex fastq files. By default, uninformative barcodes are not used to identify samples when demultiplexing fastq files (optional flag).')
  parser.add_argument('--no_expand_pcr_rows_columns', required=False, action='store_true', help='Without --no_expand_pcr_rows_columns, this program pairs the i-th P7 row with the i-th P5 column for each sample. With --no_expand_pcr_rows_columns, all combinations of P7 rows and P5 columns are used for each sample.')
  parser.add_argument('--hash_file', required=False, default=None, help='Name of hash index file. Setting hash_file to a file name enables hash index analysis.')
  parser.add_argument('-e', '--template', required=False, action='store_true', help='Write template samplesheet file (\'samplesheet.template.csv\') with standard column formats and exit (optional flag).')
  parser.add_argument('-d', '--documentation', required=False, action='store_true', help='Display documentation and exit (optional flag).')
  parser.add_argument('-v', '--version', required=False, action='store_true', help='Give program and JSON output file versions and exit (optional flag).')
  args = parser.parse_args()

  # Write documentation.
  if( args.documentation ):
    display_documentation()
    sys.exit( 0 )

  # Write versions.
  if( args.version ):
    print( 'Program version: %s' % ( program_version ) )
    print( 'JSON output file  version: %s' % ( json_file_version ) )
    sys.exit( 0 )

  # Write samplesheet template file.
  if( args.template ):
    write_samplesheet_template()
    sys.exit( 0 )

  # Check for required command line parameters.
  error_string = ''
  if( args.input == None ):
    error_string += '  input filename parameter: -i <input filename> or --input <input filename>\n'
  if( args.output == None ):
    error_string += '  output filename parameter: -o <output filename> or --output <output filename>\n'
  if( len( error_string ) > 0 ):
    print( 'Error: missing command line parameters\n%s' % ( error_string ) )
    print( 'For help run \'sciatac_samplesheet.py -h\' or \'sciatac_samplesheet.py -d\'' )
    sys.exit( -1 )

  #
  # Check command line parameter consistency.
  #
  check_args( args )

  # Go to work.
  filename_in = args.input
  filename_out = args.output

  column_name_list, samplesheet_row_list = read_samplesheet( open( filename_in, newline='' ) )

  clean_samplesheet_data( column_name_list, samplesheet_row_list )
  samplesheet_row_list = check_sample_names( column_name_list, samplesheet_row_list )
  check_genome_names( column_name_list, samplesheet_row_list )
  check_peak_groups( column_name_list, samplesheet_row_list )
  check_peak_files( column_name_list, samplesheet_row_list )
  check_external_sample_name( column_name_list, samplesheet_row_list )
  check_tissue( column_name_list, samplesheet_row_list )

  check_peak_spec( column_name_list, samplesheet_row_list )
  check_wrap_group( column_name_list, samplesheet_row_list )

  # Build a peak files dict before splitting samples. Pass the
  # peak_files_dict to write_samplesheet_json_format() and
  # write it in the output json file. At this time (20240129),
  # the peak_files map is made in the bbi-sciatac-demux main.nf
  # script and stored in demux_out/args.json. Making the 
  # peak_files_dict may make sense if we split the samples
  # by small sets of wells in the future.
#  peak_files_dict = get_peak_files(  column_name_list, samplesheet_row_list )


  if(not args.no_expand_pcr_rows_columns):
    samplesheet_row_list = expand_sample_rows(column_name_list, samplesheet_row_list)

  row_out_list = make_samplesheet_indexes( column_name_list, samplesheet_row_list )
  check_sample_identifier( row_out_list, args.sample_identifier )
  wrap_groups_dict = get_wrap_groups( row_out_list )
  if( args.format == 'json' ):
    write_samplesheet_json_format( open( filename_out, 'w' ), column_name_list, samplesheet_row_list, row_out_list, wrap_groups_dict, level = args.level, number_wells = args.number_wells, tn5_barcodes = args.tn5_barcodes, use_all_barcodes = args.use_all_barcodes, illumina_run_directory = args.run_dir, hash_file = args.hash_file )
  else:
    write_samplesheet_index_format( open( filename_out, 'w' ), row_out_list )
  samplesheet_report( samplesheet_row_list, row_out_list, wrap_groups_dict, args )
  # diagnostic dump
  # dump_row_out_list( row_out_list )

