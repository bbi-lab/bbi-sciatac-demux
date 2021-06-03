#!/usr/bin/env python3


#
# Program: pipeline_logger.py
# Purpose: write processing pipeline log file entries.
#


import os
import subprocess
import argparse
import textwrap
import json


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
  parser.add_argument('-c', '--process_command_list', nargs='+', required=False, help='List of process commands separated by whitespace.')
  parser.add_argument('-v', '--version_command_list', nargs='+', required=True, help='List of version information commands separated by whitespace.')
  parser.add_argument('-f', '--file_name_list', nargs='+', required=False, help='Names of file with additional logging information.')
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
    out_dict['additional_information'] = []
    for file_name in args.file_name_list:
      ifp = open(file_name, 'r')
      out_dict['additional_information'].append( {'file_name': file_name, 'file_content': ifp.read().rstrip()})

  # Write to JSON file.
  if(args.output_file == None):
    log_file_name = args.sample_name + '.' + args.start_time.replace(':', '') + '.' + args.process_name + '.' + args.nextflow_run_name + '.' + work_id + '.json'
  else:
    log_file_name = args.output_file
  if(args.output_directory != None):
    log_file_name = os.path.join(args.output_directory, log_file_name)

  with open(log_file_name, 'w') as json_file:
    json.dump(out_dict, json_file)

