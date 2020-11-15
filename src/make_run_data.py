#!/usr/bin/env python3


import sys
import argparse
import json
import re

def read_args_json( args_file ):
  """
  Read args.json file created in bbi-sciatac-demux pipeline run.
  """
  args_json = json.load( args_file )
  return( args_json )


def read_lane_stats( lane_file ):
  """
  Read RUN001_Lxxx.stats.json file created in bbi-sciatac-demux pipeline run.
  These files are in the demux_out/fastqs_barcode directory. The run is always
  RUN001 and the Lxxx denotes the results for each lane; that is, L001, ...
  """
  lane_stats = json.load( lane_file )
  return( lane_stats )


def make_lane_stats_list( args_json ):
  """
  Make a list of lane stats from lane stats files.
  """
  lane_stats_list = []
  num_lanes = args_json['RUN001']['sequencing_run_metadata']['lanes']
  for ilane in range( num_lanes ):
    lane_file = '%s/fastqs_barcode/RUN001_L%03d.stats.json' % ( args_json['RUN001']['DEMUX']['demux_dir'], ilane + 1 )
    lane_stats_list.append( read_lane_stats( open( lane_file ) ) )
  return( lane_stats_list )


def find_range_extrema( range_sublist, range_extrema ):
  ranges = range_sublist.split( ',' )
  min = range_extrema[0]
  max = range_extrema[1]
  for range in ranges:
    mobj = re.match( '([0-9]+)[-]([0-9]+)', range )
    if( mobj ):
      i0 = int( mobj.group( 1 ) )
      i1 = int( mobj.group( 2 ) )
      if( min is None or i0 < min ):
        min = i0
      if( max is None or i1 > max ):
        max = i1
    else:
      mobj = re.match( '([0-9]+)', range )
      if( mobj ):
        i = int( mobj.group( 1 ) )
        if( min is None or i < min ):
          min = i
        if( max is None or i > max ):
          max = i
      else:
        raise ValueError( 'bad range in sample_index_list: %s' % ( range_sublist ) )
  range_extrema[0] = min
  range_extrema[1] = max
  return( range_extrema )


def make_lig_plate_list( args_json ):
  sample_index_list = args_json['RUN001']['sample_data']['sample_index_list']
  num_sample = len( sample_index_list )
  n7_range_extrema = [ None ] * 2
  n5_range_extrema = [ None ] * 2
  for isample in range( num_sample ):
    sample_dict = sample_index_list[isample]
    ranges = sample_dict['ranges']
    range_list = ranges.split( ':' )
    n7_range_extrema = find_range_extrema( range_list[0], n7_range_extrema )
    n5_range_extrema = find_range_extrema( range_list[3], n5_range_extrema )
  range_extrema = [0] * 2
  range_extrema[0] = n7_range_extrema[0] if( n7_range_extrema[0] >= n5_range_extrema[0] ) else n5_range_extrema[0]
  range_extrema[1] = n7_range_extrema[1] if( n7_range_extrema[1] >= n5_range_extrema[1] ) else n5_range_extrema[1]
  iplate_1 = int( range_extrema[0] / 96 ) + 1
  iplate_n = int( range_extrema[1] / 96 ) + 1
  lig_plate_list = []
  for iplate in range( iplate_1, iplate_n ):
    lig_plate_list.append( 'P%d' % ( iplate ) )
  return( lig_plate_list )


def make_pcr_combo_list( args_json ):
  if( args_json['RUN001']['sample_data']['pcr_format' ] != 'row_col' ):
    return( None )
  p7_row_list = args_json['RUN001']['sample_data']['p7_row_list']
  p5_col_list = args_json['RUN001']['sample_data']['p5_col_list']
  assert len( p7_row_list ) == len( p5_col_list ), 'different lengths for p7_row_list and p5_col_list'
  pcr_combo_list = []
  for ( pcr_p7, pcr_p5 ) in zip( p7_row_list, p5_col_list ):
    pcr_combo_list.append( '%s%s' % ( pcr_p7.upper(), pcr_p5.upper() ) )
  return( pcr_combo_list )


def make_run_data_dict( args_json ):
  run_name_list = ['sci-ATAC-seq']
  lane_stats_list = make_lane_stats_list( args_json )
  num_lanes = args_json['RUN001']['sequencing_run_metadata']['lanes']
  lane_list = []
  for ilane in range( num_lanes ):
    lane_list.append( '%d' % ( ilane + 1 ) )
  level_list = []
  level_list.append( '%d' % ( args_json['RUN001']['DEMUX']['level'] ) )
  lig_plate_list = make_lig_plate_list( args_json )
  pcr_combo_list = make_pcr_combo_list( args_json )
  lig_combo_list = lig_plate_list
  bad_well_barcodes = {}
  bad_wells = []
  run_data_dict = { 'run_name' : run_name_list,
                    'lane_list' : lane_list,
                    'plate_list' : lig_plate_list,
                    'pcr_combo_list' : pcr_combo_list,
                    'lig_combo_list' : lig_combo_list,
                    'level' : level_list,
                    'lane_stats' : lane_stats_list,
                    'include_norm' : True,
                    'bad_wells_barcodes' : bad_well_barcodes,
                    'bad_wells' : bad_wells
                  }
  return( run_data_dict )


def write_run_data( run_data_dict, out_file ):
  json_object = json.dumps( run_data_dict, indent = 2 )
  out_file.write( 'const run_data =\n' )
  out_file.write( json_object )
  return( 0 )


if __name__ == '__main__':

  parser = argparse.ArgumentParser(description='A program to write run_data.js file for sci-ATAC demux dashboard.')
  parser.add_argument('-i', '--input_file', required=True, help='Path to args.json input file in sci-ATAC demux output directory.')
  parser.add_argument('-o', '--output_file', required=True, help='Name of output file.')
  args = parser.parse_args()

  args_json = read_args_json( open( args.input_file ) )
  run_data_dict = make_run_data_dict( args_json )
  write_run_data( run_data_dict, open( args.output_file, 'w' ) )


