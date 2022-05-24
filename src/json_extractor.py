#!/usr/bin/env python3

# Extract blocks from JSON file.

# The official repository for the program is
#
#   https://github.com/bbi-lab/bbi-sciatac-demux/src/json_extractor.py
#
# The program is distributed from there to
#
#   https://github.com/bbi-lab/bbi-sciatac-analyze/src/json_extractor.py
#

import os
import sys
import subprocess
import argparse
import textwrap
import json


json_extractor_version = '20220429a'


# Read a process-level log file and return the contents
# as a dictionary. The file contents are in JSON format
# so this function simply reads the file into a dictionary.
def read_json_file(file_name):
  json_struct = json.load(open(file_name))
  return(json_struct)


def merge_dicts(dict_src, dict_dst):
  for key_src in dict_src.keys():
    if(not key_src in dict_dst.keys()):
      dict_dst[key_src] = dict_src[key_src]
  return(dict_dst)


def json_extract(json_obj, key_list):
    json_blocks = {}
    for key in json_obj.keys():
      if(key in key_list):
        json_blocks[key] = json_obj[key]
        continue
      if(type(json_obj[key]) == dict):
        json_blocks_sub = json_extract(json_obj[key], key_list)
        json_blocks = merge_dicts( json_blocks_sub, json_blocks)
    return json_blocks


# Indent text block but not the leading line.
def indent_text_block(text_block, spacer):
  parts = text_block.split('\n')
  indent = '\n' + spacer
  out_text_block = indent.join(parts)
  return(out_text_block)
  

# Write a list of data. The list may consist
# of strings, lists, and dictionaries.
def write_data_list(data_list, indent_level=0, output_file=sys.stdout):
  spacer = ' ' * indent_level * 2
  for item in data_list:
    if(type(item) is dict):
      print('%s+ dictionary item'% (spacer), file=output_file)
      write_data_dict(item, indent_level+1, output_file)
    elif(type(item) is list):
      print('%s+ list item'% (spacer), file=output_file)
      write_data_list(item, indent_level+1, output_file)
    else:
      if(type(item) is str):
        out_text_block = indent_text_block(item.rstrip(), ' ' * 2 * (indent_level+1))
      else:
        out_text_block = str(item)
      print('%slist item: %s' % (spacer, out_text_block), file=output_file)


# Write a dictionary of data. The dictionary
# may consist of strings, lists, and
# dictionaries.
def write_data_dict(data_dict, indent_level=0, output_file=sys.stdout):
  for key in data_dict.keys():
    spacer = ' ' * indent_level * 2
    if(type(data_dict[key]) is dict):
      print('%s+ %s:' % ( spacer, key), file=output_file)
      write_data_dict(data_dict[key], indent_level+1, output_file)
    elif(type(data_dict[key]) is list):
      print('%s+ %s:' % ( spacer, key), file=output_file)
      write_data_list(data_dict[key], indent_level+1, output_file)
    else:
      if(type(data_dict[key]) is str):
        out_text_block = data_dict[key].rstrip()
        if('\n' in out_text_block):
          out_text_block = '\n' + out_text_block
        out_text_block = indent_text_block(out_text_block, ' ' * 2 * (indent_level+1))
      else:
        out_text_block = str(data_dict[key])
      print('%s+ %s: %s' % (spacer, key, out_text_block), file=output_file)


#json_struct = read_json_file(file_name)

#json_sub = json_extract(json_struct, ['ReadInfosForLanes', 'ConversionResults'])

# json_struct['ConversionResults'])

#write_data_dict(json_sub, indent_level=0)


if __name__ == '__main__':
  parser = argparse.ArgumentParser('Script to extract dictionary blocks from JSON file.')
  parser.add_argument('-i', '--input_files', required=True, nargs='+', default=None, help='A list of input JSON files.')
  parser.add_argument('-k', '--keys', required=True, nargs='+', default=None, help='A list of keys of blocks to extract.')
  parser.add_argument('-o', '--output_file', required=False, default=None, help='Optional JSON output filename. Default to stdout.')
  args = parser.parse_args()

  out_dict ={}
  for file_name in args.input_files:
    json_obj = json_extract(read_json_file(file_name), args.keys)
    out_dict['File name: ' + file_name] = json_obj

  if(not args.output_file is None):
      with open(args.output_file, 'w') as json_file:
        json.dump(out_dict, json_file)
  else:
    print(out_dict)



