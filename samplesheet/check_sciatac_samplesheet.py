#!/usr/bin/env python3

#
# Compare json well indices to input samplesheet well specs.
#

import sys
import re
import csv
import json
import argparse


def make_column_index(input_samplesheet_column_names):
  column_index_dict = {}
  for i, column_dict in enumerate(input_samplesheet_column_names):
    column_index_dict[column_dict['type']] = [i, column_dict['format']]
  return(column_index_dict)


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

  zero_pad_col = True
  id_length = 2
  well_id = 'P%d-%s%s' % (ipl + 1, well_row, pad_well_col(well_col, zero_pad_col, id_length))
#  well_id = '%s%s' % ( well_row, pad_well_col( well_col, True, 2 ) )

  return( well_id )


def convert_index_string(in_index, across_row_first, format):
  if(format == 'indexes'):
    return(in_index)

  # Split index groups.
  out_string = ''
  index_groups_list = in_index.split(',')
  for igroup, index_group in enumerate(index_groups_list):
    if(igroup > 0):
      out_string = out_string + ','
    # Split ranges.
    index_list = index_group.split('-')
    if(len(index_list) == 1):
      well_id = index_to_well(int(index_list[0])-1, across_row_first)
      out_string = out_string + well_id
    elif(len(index_list) == 2):
      well_id = index_to_well(int(index_list[0])-1, across_row_first)
      out_string = out_string + well_id + ':'
      well_id = index_to_well(int(index_list[1])-1, across_row_first)
      out_string = out_string + well_id
    else:
      print('Error: unexpected number of tokens in range')
      sys.exit(-1)
  return(out_string)


if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='A program to check sci-ATAC json samplesheet file.')

  parser.add_argument('-i', '--input', required=True, default=None, help='Input json samplesheet filename (required string).')
  args = parser.parse_args()

  fp = open(args.input)
  json_data = json.load(fp)

  print('input column name dicts (len: %d): ' % (len(json_data['input_samplesheet_column_names'])), json_data['input_samplesheet_column_names'])

  column_index_dict = make_column_index(json_data['input_samplesheet_column_names'])

  print('%s %d  %s %d  %s %d  %s %d' % ('n5', column_index_dict['n5'][0],
                                        'n7', column_index_dict['n7'][0],
                                        'p5', column_index_dict['p5'][0],
                                        'p7', column_index_dict['p7'][0]))

  number_rows = len(json_data['input_samplesheet_rows'])
  print('number rows: %d' % (number_rows))

  for irow in range(number_rows):
    input_row = json_data['input_samplesheet_rows'][irow]
    print('input row: ', input_row)
    cols = list(csv.reader(input_row, delimiter=','))
    toks = []
    for col in cols:
      if(len(col) == 2 and col[0] == '' and col[1] == ''):
        continue
      toks.append(col)
#    print('cols: ', cols)
#    print('row toks: ', toks)
    row_index_list = json_data['sample_index_list'][irow]['ranges'].split(':')
    print('sample: %s' % (json_data['sample_index_list'][irow]['sample_id']))
    print('n7:')
    print('  csv:  %s' % (toks[column_index_dict['n7'][0]][0]))
    print('  json: %s' % (convert_index_string(row_index_list[0], True, column_index_dict['n7'][1])))
    print('  json: %s' % (row_index_list[0]))
    print('p7:')
    print('  csv:  %s' % (toks[column_index_dict['p7'][0]][0]))
    print('  json: %s' % (convert_index_string(row_index_list[1], True, column_index_dict['p7'][1])))
    print('  json: %s' % (row_index_list[1]))
    print('p5:')
    print('  csv:  %s' % (toks[column_index_dict['p5'][0]][0]))
    print('  json: %s' % (convert_index_string(row_index_list[2], False, column_index_dict['p5'][1])))
    print('  json: %s' % (row_index_list[2]))
    print('n5:')
    print('  csv:  %s' % (toks[column_index_dict['n5'][0]][0]))
    print('  json: %s' % (convert_index_string(row_index_list[3], False, column_index_dict['n5'][1])))
    print('  json: %s' % (row_index_list[3]))
    print()

