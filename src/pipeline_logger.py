#!/usr/bin/env python3


#
# Program: pipeline_logger.py
# Purpose: write processing pipeline log file entries.
#

# Notes:
#   o  getting JSON strings from the Nextflow process block to
#      this pipeline_logger.py program is a challenge. I managed
#      to get it to work using
#
#        MITO_READS="{\\\"mitochondrial_reads\\\": ${alignMap['has_whitelist_with_mt']}}"
#        pipeline_logger.py ... -j "\${MITO_READS}" ...
#     
#      Maybe there is an easier way?
#
# The file names appear to be complex but the complexity ensures that
# the files are not over-written by different pipeline runs, different
# process blocks, different samples, and, in some process blocks, more
# than one block run for the same sample, process block, and pipeline
# run. An example, is the callMotifsProcess in the bbi-sciatac-analyze
# pipeline where the block is executed for each gc_bin (for a sample
# and the process block). The run work directory hash in the file
# name ensures that the log file names are distinct.
#
# The official repository for the program is
#
#   https://github.com/bbi-lab/bbi-sciatac-demux/src/pipeline_logger.py
#
# The program is distributed from there to
#
#   https://github.com/bbi-lab/bbi-sciatac-analyze/src/pipeline_logger.py
#

import os
import subprocess
import argparse
import textwrap
import json

pipeline_logger_version = '20220429c'

def split_path(path):
  split_parts = []
  while 1:
    parts = os.path.split(path)
    if parts[0] == path:
      split_parts.insert(0, parts[0])
      break
    elif parts[1] == path:
      split_parts.insert(0, parts[1])
      break
    else:
      path = parts[0]
      split_parts.insert(0, parts[1])
  return split_parts


if __name__ == '__main__':
  parser = argparse.ArgumentParser('Script to log processing pipeline run information.')
  parser.add_argument('-r', '--nextflow_run_name', required=True, help='Nextflow run name referenced by \'$workflow.runName\'.')
  parser.add_argument('-n', '--sample_name', required=True, help='Sample name.')
  parser.add_argument('-p', '--process_name', required=True, help='Nextflow process name (prefix the name with a block sequence number in order to sort files).')
  parser.add_argument('-x', '--extended_identifier', required=False, default=None, help='Additional identifier used when a sample + process combination occurs more than once per run.')
  parser.add_argument('-c', '--process_command_list', required=False, nargs='*', help='List of process commands separated by whitespace.')
  parser.add_argument('-v', '--version_command_list', required=False, nargs='*', help='List of version information commands separated by whitespace.')
  parser.add_argument('-f', '--file_name_list', required=False, nargs='*', help='Names of file with additional logging information.')
  parser.add_argument('-t', '--file_name_list_trimmed', required=False, nargs='*', help='Names of file with additional logging information where the first token is to be trimmed off.')
  parser.add_argument('-j', '--json_strings', required=False, nargs='*', help='JSON string.')
  parser.add_argument('-J', '--json_files', required=False, nargs='*', help='JSON files.')
  parser.add_argument('-s', '--start_time', required=True, help='Process start time.')
  parser.add_argument('-e', '--end_time', required=True, help='Process end time.')
  parser.add_argument('-d', '--output_directory', required=False, help='Output directory.')
  parser.add_argument('-o', '--output_file', required=False, help='Optional output filename.')
  args = parser.parse_args()

  # Directory in which this script runs.
  cwd = os.getcwd()
  cwd_parts = split_path(cwd)
  cwd_parts[-1] = cwd_parts[-1][0:6]
  nparts = len( cwd_parts )
  work_id = '_'.join(cwd_parts[-2:])

  # Save to dictionary in order to write to json file.
  out_dict = {}

  # Process run information.
  out_dict['process'] = {}
  out_dict['process']['run_name'] = args.nextflow_run_name
  out_dict['process']['sample_name'] = args.sample_name
  if(not args.extended_identifier is None):
    out_dict['process']['extended_identifier'] = args.extended_identifier
  out_dict['process']['process_name'] = args.process_name
  out_dict['process']['work_directory'] = cwd
  out_dict['process']['start_time'] = args.start_time
  out_dict['process']['end_time'] = args.end_time


  # Run version commands and collect output.
  out_dict['program_versions'] = []
  for command in args.version_command_list:
    out_dict['program_versions'].append({'command': command.rstrip(),
                                         'output': subprocess.run(command, shell=True, stderr=subprocess.STDOUT, stdout=subprocess.PIPE).stdout.decode('ascii')})

  # Process commands.
  if(args.process_command_list != None):
    out_dict['command_list'] = []
    for command in args.process_command_list:
      out_dict['command_list'].append(command.rstrip())

  # Additional logged text.
  if(args.file_name_list != None):
    out_dict['file_contents'] = []
    for file_name in args.file_name_list:
      if(os.path.exists(file_name)):
        ifp = open(file_name, 'r')
        out_dict['file_contents'].append( {'file_name': file_name, 'file_content': ifp.read().rstrip()} )

  # Read file contents, trim off first token, and
  # save a json dictionary entry.
  # This is to remove the TAGs on the bbi-genome-data
  # record logs.
  if(args.file_name_list_trimmed != None):
    out_dict['file_contents_trimmed'] = []
    for file_name in args.file_name_list_trimmed:
      if(os.path.exists(file_name)):
        ifp = open(file_name, 'r')
        in_lines = ifp.readlines()
        for i, in_line in enumerate(in_lines):
          in_tokens = in_line.split(' ')
          in_lines[i] = ' '.join(in_tokens[1:])
        in_text = ''.join(in_lines)
        out_dict['file_contents_trimmed'].append( {'file_name': file_name, 'file_content_trimmed': in_text} )

  # Store JSON strings.
  if(args.json_strings != None):
    out_dict['json_strings'] = []
    for json_string in args.json_strings:
      out_dict['json_strings'].append(json.loads(json_string))

  # Store the contents of JSON files.
  if(args.json_files != None):
    out_dict['json_files'] = []
    for json_file in args.json_files:
      json_obj = json.load(open(json_file))
      out_dict['json_files'].append({'file_name': json_file, 'file_content': json_obj})

  # Write to JSON file.
  if(args.output_file == None):
    log_file_name = args.sample_name + '.' + args.start_time.replace(':', '') + '.' + args.process_name + '.' + args.nextflow_run_name + '.' + work_id + '.json'
  else:
    log_file_name = args.output_file
  if(args.output_directory != None):
    log_file_name = os.path.join(args.output_directory, log_file_name)

  with open(log_file_name, 'w') as json_file:
    json.dump(out_dict, json_file)

