#!/usr/bin/python3

#
# Make a sample sheet suitable for Andrew Hill's sci-ATAC-seq pipeline from
# a BBI sample sheet file that has the format
#
# <ligation barcode well id>,<sample_name>,<genome_name>
#
# Note: this description is now incomplete...
#


import re
import sys
import csv
import json
import argparse


#
# This genomes map will require additions, modifications, or, possibly, removal.
#
genomes_map = {
 'Human': 'human',
 'Mouse': 'mouse',
 'Barnyard': 'barnyard'
}


#
# Notes:
#   o  this script converts well to the indices used by Andrew's
#      (downstream) scripts.
#   o  the README.md file in Andrew's sciatac-pipeline repository
#      has the following description of some parameters and the
#      expected samplesheet file
#      # Demux
#      The first step is to demux your sequencing run(s).
#      
#      ```
#      ./demux_sciatac \
#              --rundir /net/shendure/vol9/seq/NEXTSEQ/180823_NS500488_0677_AH3WGGAFXY \
#              --outdir demux_run1 \
#              --samplesheet sample_sheet.txt \
#              --chemistry 3level_384
#      ```
#      
#      - `--rundir` is the BCL directory for the sequencing run
#      - `--outdir` is the output directory for this run of the pipeline (one per seq run)
#      - `--nthreads` sets the number of threads
#      - `--samplesheet` is a tab-delimited file with the following format:
#      ```
#      sample_id      ranges
#      sample_1       1-96:25-48:65-80:1-20
#      sample_2       1-96:25-48:65-80:21-40
#      sample_3       1-96:25-48:65-80:41-60
#      sample_4       1-96:25-48:65-80:61-80
#      sample_5       1-96:25-48:65-80:80-96
#      ```
#      
#      The `sample_id` column can have whatever you want as an ID for each sample in the run.
#      
#      The `ranges` column specifies the range of `<N7 ligation indices>:<N7 PCR indices>:<N5 PCR indices>:<N5 ligation indices>` used for this run. Note that the N7/P7 indices are numbered by column and the N5/P5 indices are numbered by row. This format is different for two-level (barcoded Tn5) chemistry.
#      
#      This is used to make a table of all the possible indices for each sample, and then used in the remainder of the pipeline.
#      
#      - `--chemistry` default assay is `3level_384`, which corresponds to ligation based assay with plates from set of 384 (such as human fetal atlas). `3level_96` is also available if using the original set of 96-well barcode designs (most of 3level ligation-based assay development). `2level` is also available for the original barcoded Tn5-based assay. 
#      
#      You may repeat this for each sequencing run that you have for a given set of samples.

def well_str_to_list( wells ):
  """
  Convert a string of well tokens (e.g. "A05 B14 G08") to a list
  (0-pad column token).

  Args:
    wells: string of well tokes (e.g. "A05 B14 G08")

  Returns:
    well_str_list: a list of well tokens

  """
  well_str_list = []
  for well in wells.split():
    mobj = re.match( '^([A-H])([0-9]{1,2}$)', well )
    if( mobj == None ):
      print( 'Error: bad well name %s' % ( well_str ), file=sys.stderr )
      sys.exit( -1 )
    well_str_list.append( '%s%02d' % ( mobj.group(1), int( mobj.group(2) ) ) )
  return( well_str_list )


def well_to_index( well_str, row_ordered = True ):
  """
  Convert a well string to a plate index (well A01=1).

  Args:
    well_str: a string of a well (e.g. A10)
    row_ordered: True=increases by row, i.e., A01=1, B01=2, ... (n5/p5 are row_ordered)

  Returns:
    index: an integer well index (1-96)

  """
  mobj = re.match( '^([A-H])([0-9]{1,2})$', well_str )
  if( mobj == None ):
    print( 'Error: bad well string \'%s\'' % ( well_str ), file=sys.stderr )
    sys.exit( -1 )
  nrow = ord( mobj.group( 1 ) ) - ord( 'A' )
  ncol = int( mobj.group( 2 ) ) - 1
  if( row_ordered ):
    well_index = ncol * 8 + nrow + 1
  else:
    well_index = nrow * 12 + ncol + 1
  return( well_index )


def plate_well_to_index( plate_well_str, row_ordered = True ):
  """
  Convert a plate-well string to a plate index (well P1-A01=1).

  Args:
    plate_well_str: a string of a plate-well (e.g. P1-A10)
    row_ordered: True=increases by row, i.e., A01=1, B01=2, ... (n5/p5 are row_ordered)

  Returns:
    index: an integer well index (1-96)

  """
  mobj = re.match( '^P([01]?[1-4])-([A-H][0-9]{1,2})$', plate_well_str )
  if( mobj == None ):
    print( 'Error: bad plate-well string \'%s\'' % ( plate_well_str ), file=sys.stderr )
    sys.exit( -1 )
  nplate = int( mobj.group( 1 ) ) - 1
  well_str = mobj.group( 2 )
  plate_well_index = nplate * 96 + well_to_index( well_str, row_ordered )
  return( plate_well_index )


def read_samplesheet( samplesheet_file ):
  """
  Read and parse a BBI sample sheet CSV file.

  Args:
    samplesheet_file: a string with the path to the samplesheet CSV file

  Returns:
    samples (set)
    genomes_dict  dictionary of genome names keyed by sample name
    well_samples  (list) well sample names
    well_genomes  (list) well genomes
    well_n7_indices (list) N7 well indices
    well_n5_indices (list) N5 well indices

  """
  well_samples = []
  well_genomes = []
  well_n7_indices = []
  well_n5_indices = []
  row_list = []
 
  with open( samplesheet_file ) as fp:
    csv_reader = csv.reader( fp, delimiter=',' )
    for row in csv_reader:
      row_list.append(row)
      if( len( row[0] ) == 0 ):
        break
#      print( '%s -> %s    %s -> %s  %s' % ( row[0], n7_well_to_index[row[0]], row[0], n5_well_to_index[row[0]], row[1] ) )
      well_samples.append( row[1] )
      well_genomes.append( row[2] )
      well_n7_indices.append( plate_well_to_index( row[0], row_ordered = False ) )
      well_n5_indices.append( plate_well_to_index( row[0], row_ordered = True ) )
  
  samples = set( well_samples )
  
  genomes_dict = {}
  for i in range( len( well_samples ) ):
    genomes_dict.setdefault( well_samples[i], well_genomes[i] )

  return( samples, genomes_dict, well_samples, well_genomes, well_n7_indices, well_n5_indices, row_list )


def make_tag_masks( samples, well_samples, well_n7_indices, well_n5_indices ):
  """
  Make masks of ligation sequence wells used in run.
  return
    n7_masks  (list)
    n5_masks  (list )
  """
  n7_masks = {}
  for sample in samples:
    n7_masks[sample] = [0] * 384
  
  n5_masks = {}
  for sample in samples:
    n5_masks[sample] = [0] * 384
  
  
  nwells = len( well_samples )
  for iwell in range( nwells ):
    sample = well_samples[iwell]
    assert n7_masks[sample][well_n7_indices[iwell]-1] == 0, 'duplicate sample wells in sample sheet'
    n7_masks[sample][well_n7_indices[iwell]-1] = 1
    assert n5_masks[sample][well_n5_indices[iwell]-1] == 0, 'duplicate sample wells in sample sheet'
    n5_masks[sample][well_n5_indices[iwell]-1] = 1

  return( n7_masks, n5_masks )


def make_pcr_ranges_rows_cols( pcr_rows, pcr_cols ):
  """
  Make a string of i7 and i5 PCR well indices used from rows+cols.
  return
    pcr_ranges_str  (string )
  """
  beg = 0
  end = 0
  i7_mask = [0] * 96
  i5_mask = [0] * 96
  for row in pcr_rows:
    row = row.upper()
    beg = ( ord( row ) - ord( 'A' ) ) * 12
    end = beg + 11
    for i in range( beg, end + 1 ):
      i7_mask[i] = 1
  for col in pcr_cols:
    beg = ( int( col ) - 1 ) * 8
    end = beg + 7
    for i in range( beg, end + 1 ):
      i5_mask[i] = 1

  i7_ranges = []
  flag = 0
  for iwell in range( 96 ):
    if( i7_mask[iwell] and flag == 0 ):
      beg = iwell
      flag = 1
    elif( ( not i7_mask[iwell] ) and flag ):
      end = iwell - 1 
      if( beg != end ):
        i7_ranges.append( '%d-%d' % ( beg + 1, end + 1) )
      else:
        i7_ranges.append( '%d' % ( beg + 1) )
      flag = 0
  if( flag ):
    if( ( beg + 1 ) != 96 ):
      i7_ranges.append( '%d-%d' % ( beg + 1, 96 ) )
    else:
      i7_ranges.append( '%d' % ( 96 ) )

  i5_ranges = []
  flag = 0
  for iwell in range( 96 ):
    if( i5_mask[iwell] and flag == 0 ):
      beg = iwell
      flag = 1
    elif( ( not i5_mask[iwell] ) and flag ):
      end = iwell - 1
      if( beg != end ):
        i5_ranges.append( '%d-%d' % (beg + 1, end + 1) )
      else:
        i5_ranges.append( '%d' % (beg + 1) )
      flag = 0
  if( flag ):
    if( ( beg + 1 ) != 96 ):  
      i5_ranges.append( '%d-%d' % ( beg + 1, 96 ) )
    else:
      i5_ranges.append( '%d' % (96) )

  pcr_ranges_str = ''
  for i, index_range in enumerate( i7_ranges ):
    if( i > 0 ):
      pcr_ranges_str += ','
    pcr_ranges_str += index_range

  pcr_ranges_str += ':'

  for i, index_range in enumerate( i5_ranges ):
    if( i > 0 ):
      pcr_ranges_str += ','
    pcr_ranges_str += index_range

  return( pcr_ranges_str )


def range_extract( lst ):
  """
  Convert a list of integers to tuples of integer ranges.

  Args:
    lst: list of integers

  Returns:
    Yield 2-tuple ranges or 1-tuple single elements from list of increasing ints'

  """
  lst.sort()
  lenlst = len( lst )
  i = 0
  while( i < lenlst ):
    low = lst[i]
    while( i < lenlst - 1 and lst[i] + 1 == lst[i+1] ):
      i += 1
    hi = lst[i]
    if( hi - low >= 1 ):
      yield( low, hi )
    else:
      yield( low, )
    i += 1
 

def ranges_to_string( ranges ):
  """
  Convert range tuples to string

  Args:
    ranges: tuples of integer ranges

  Returns:
    string of ranges

  """
  return( ','.join( ( ( '%i-%i' % r ) if len( r ) == 2 else '%i' % r ) for r in ranges ) )
 

def make_pcr_ranges_wells( p7_wells, p5_wells ):
  """
  Convert PCR well strings to index string.

  Args:
    p7_wells: list of p7 well strings
    p5_wells: list of p5 well strings

  Returns:
    pcr_ranges_str: string of PCR range indices

  """
  p7_index = []
  p5_index = []
  for well_str in p7_wells:
    p7_index.append( well_to_index( well_str, row_ordered = False ) )
  for well_str in p5_wells:
    p5_index.append( well_to_index( well_str, row_ordered = True ) )
  p7_index = list( set( p7_index ) )
  p5_index = list( set( p5_index ) )
  pcr_ranges_str = ranges_to_string( range_extract( p7_index ) ) + ':' + ranges_to_string( range_extract( p5_index ) )
  return( pcr_ranges_str )


def check_pcr_parameters( args ):
  """
  Check if PCR wells are given by row+col (well_wise=False) or
  by wells (well_wise=True).

  Args:
    args: argsparse

  Returns:
    well_wise: boolean True=well-wise PCR wells; False=row+col PCR wells
    error_flag: integer 0=OK; -1=Error

  """
  error_flag = 0
  if( args.p7_rows != 0 and args.p7_rows != 'none' and
      args.p5_cols != 0 and args.p5_cols != 'none' and
      ( args.p7_wells == 0 or args.p7_wells == 'none' and
        args.p5_wells == 0 or args.p5_wells == 'none' ) ):
     well_wise = False
     # by row+col...
     # check rows
     str_list = args.p7_rows.split()
     nrow = len( str_list )
     for str in str_list:
       if( re.search( '^[A-H]$', str ) == None ):
         print( 'Error: bad p7 row list' )
         error_flag = -1
         break
     # check cols
     str_list = args.p5_cols.split()
     ncol = len( str_list )
     for str in str_list:
       if( re.search( '^[0-9]{1,2}$', str ) == None ):
         print( 'Error: bad p5 col list' )
         error_flag = -1
         break
     # check that nrow == ncol
     if( nrow != ncol ):
       print( 'Error: inconsistent number of PCR rows (%d) and cols (%d)' % ( nrow, ncol ) )
       error_flag = -1
  elif( ( args.p7_rows == 0 or args.p7_rows == 'none' and
           args.p5_cols == 0 or args.p5_cols == 'none' ) and
         args.p7_wells != 0 and args.p7_wells != 'none' and
         args.p5_wells != 0 and args.p5_wells != 'none' ):
    # by wells
    well_wise = True
    well_list = well_str_to_list( args.p7_wells )
    np7_well = len( well_list )
    well_list = well_str_to_list( args.p5_wells )
    np5_well = len( well_list )
    if( np7_well != np5_well ):
       print( 'Error: inconsistent number of PCR p7 wells (%d) and p5 wells (%d)' % ( np7_well, np5_well ) )
       error_flag = -1
  else:
    print( 'Error: mixed PCR parameters: p7_rows and p5_rows and p7_wells and p5_wells' )
    error_flag = -1
  if( error_flag == -1 ):
    well_wise = False
  return( well_wise, error_flag )


def make_pcr_ranges( args ):
  """
  Make PCR index string from args.

  Args:
    args: argsparse

  Returns:
    pcr_well_wise: boolean False=row+col; True=wells
    pcr_p7_list: string of p7 rows or wells
    pcr_p5_list: string of p5 cols or wells
    pcr_ranges_str: string of PCR indexes
    
  """
  ( well_wise, error_flag ) = check_pcr_parameters( args )
  if( error_flag ):
    print( 'Error: bad PCR parameters' )
    sys.exit( -1 )
  if( not well_wise ):
    p7_rows = args.p7_rows.split()
    p5_cols = args.p5_cols.split()
    pcr_ranges_str = make_pcr_ranges_rows_cols( p7_rows, p5_cols )
    return( well_wise, p7_rows, p5_cols, pcr_ranges_str )
  else:
    p7_wells = well_str_to_list( args.p7_wells )
    p5_wells = well_str_to_list( args.p5_wells )
    pcr_ranges_str = make_pcr_ranges_wells( p7_wells, p5_wells )
    return( well_wise, p7_wells, p5_wells, pcr_ranges_str )


if __name__ == '__main__':

  parser = argparse.ArgumentParser(description='A program to prepare sample index data for a run.')
  parser.add_argument('-l', '--level', type=str, default=None, required=True, help='Experiment level. Either 2 or 3.')
  parser.add_argument('-r', '--p7_rows', type=str, default=0, required=False, help='PCR plate p7 rows used for indexing entered as a single string')
  parser.add_argument('-c', '--p5_cols', type=str, default=0, required=False, help='PCR plate p5 columns used for indexing entered as a single string.')
  parser.add_argument('-7', '--p7_wells', type=str, default=0, required=False, help='PCR plate p7 wells used for indexing entered as a single string.')
  parser.add_argument('-5', '--p5_wells', type=str, default=0, required=False, help='PCR plate p5 wells used for indexing entered as a single string.')
  parser.add_argument('-s', '--samplesheet', type=str, default=None, required=True, help='Sample sheet file name (CSV file).')
  parser.add_argument('-o', '--outfile', type=str, default=None, required=True, help='Output file name.' )
  args = parser.parse_args()

  level = args.level
  samplesheet_file = args.samplesheet
  (samples, genomes_dict, well_samples, well_genomes, well_n7_indices, well_n5_indices, row_list) = read_samplesheet(samplesheet_file)
  
  (n7_masks, n5_masks) = make_tag_masks(samples, well_samples, well_n7_indices, well_n5_indices)
 
  (well_wise, error_flag) = check_pcr_parameters(args)
  if( error_flag ):
    sys.exit( -1 )

  ( pcr_well_wise, pcr_p7_list, pcr_p5_list, pcr_ranges_str ) = make_pcr_ranges( args )
  
  # print( 'N7' )
  # for sample in samples:
  #   print( 'sample: %s' % ( sample ) )
  #   iwell = 0
  #   for i in range( 16 ):
  #     for j in range( 24 ):
  #       print( ' %d' % ( n7_masks[sample][iwell] ), end='' )
  #       iwell += 1
  #     print()
  #   print()
  
  
  # print( 'N5' )
  # for sample in samples:
  #   print( 'sample: %s' % ( sample ) )
  #   iwell = 0
  #   for i in range( 16 ):
  #     for j in range( 24 ):
  #       print( ' %d' % ( n5_masks[sample][iwell] ), end='' )
  #       iwell += 1
  #     print()
  #   print()
  
 
  sample_indices_list = [] 
#  print( 'sample_id\tranges\tgenome' )
  for sample in samples:
    n7_beg = 0
    n7_end = 0
    n5_beg = 0
    n5_end = 0
    n7_ranges = []
    n5_ranges = []
    n7_flag = 0
    n5_flag = 0
    for iwell in range(384):
      if(n7_masks[sample][iwell] and n7_flag == 0):
        n7_beg = iwell + 1
        n7_flag = 1
      elif(not n7_masks[sample][iwell] and n7_flag == 1):
        n7_end = iwell
        if(n7_beg != n7_end):
          n7_ranges.append('%d-%d' % (n7_beg, n7_end))
        else:
          n7_ranges.append('%d' % (n7_beg))
  
        n7_flag = 0
    if(n7_flag):
      if(n7_beg != 384):
        n7_ranges.append('%d-%d' % (n7_beg, 384))
      else:
        n7_ranges.append('%d' % (384))
  
    for iwell in range(384):
      if(n5_masks[sample][iwell] and n5_flag == 0):
        n5_beg = iwell + 1
        n5_flag = 1
      elif(not n5_masks[sample][iwell] and n5_flag == 1):
        n5_end = iwell
        if(n5_beg != n5_end):
          n5_ranges.append('%d-%d' % (n5_beg, n5_end))
        else:
          n5_ranges.append('%d' % (n5_beg))
        n5_flag = 0
    if(n5_flag):
      if(n5_beg != 384):
        n5_ranges.append('%d-%d' % (n5_beg, 384))
      else:
        n5_ranges.append('%d' % (384))
  
  #  print('sample: %s  n7 ranges:' % (sample))
  #  for arange in n7_ranges:
  #    print( '  %s' % ( arange ) )
  #  print( 'sample: %s  n5 ranges:' % ( sample ) )
  #  for arange in n5_ranges:
  #    print( '  %s' % ( arange ) )
    n7_ranges_str = ''
    for i in range(len( n7_ranges)):
      if(i > 0):
        n7_ranges_str += ','
      n7_ranges_str += n7_ranges[i]
  #  print( 'n7 ranges: %s' % ( n7_ranges_str ) )
  
    n5_ranges_str = ''
    for i in range(len(n5_ranges)):
      if(i > 0):
        n5_ranges_str += ','
      n5_ranges_str += n5_ranges[i]
  #  print( 'n5 ranges: %s' % ( n5_ranges_str ) )
  
#    print( '%s\t%s:%s:%s\t%s' % ( sample, n7_ranges_str, pcr_ranges_str, n5_ranges_str, genomes_map[genomes_dict[sample]] ) ) 
    sample_indices_list.append( {'sample_id' : sample,
                               'ranges' : ':'.join([n7_ranges_str, pcr_ranges_str, n5_ranges_str]),
                               'genome' : genomes_map[genomes_dict[sample]] })

    sample_data = { 'level': level,
                    'samplesheet' : row_list,
                    'pcr_well_wise': pcr_well_wise,
                    'pcr_p7_list' : pcr_p7_list,
                    'pcr_p5_list' : pcr_p5_list,
                    'sample_indices_list' : sample_indices_list
                  }

    with open(args.outfile, 'wt') as f:
      f.write(json.dumps(sample_data, indent=4))


