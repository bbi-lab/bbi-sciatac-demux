#!/usr/bin/env python3


#
# Program: pipeline_logger.py
# Purpose: write processing pipeline log file entries.
#


import os
import subprocess
import argparse
import textwrap


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
  parser.add_argument('-f', '--file_name', required=False, help='Name of file with additional logging information.')
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

  # Run commands and collect output.
  out_text = []
  for command in args.version_command_list:
#    out_text.append(os.popen(command).read().rstrip())
    out_text.append(subprocess.run(command, shell=True, stderr=subprocess.STDOUT, stdout=subprocess.PIPE).stdout.decode('ascii'))

  log_text = ''
  log_text += '========================================\n'
  log_text += '\n'
  log_text += 'Run name: %s\n' % (args.nextflow_run_name)
  log_text += 'Sample: %s\n' % (args.sample_name)
  log_text += 'Process: %s\n' % (args.process_name)
  log_text += 'Work directory: %s\n' % (cwd)
  log_text += 'Start time: %s\n' % (args.start_time)
  log_text += 'End time: %s\n' % (args.end_time)

  if(args.process_command_list != None):
    log_text += 'Process commands:\n'
    for command in args.process_command_list:
      log_text += '  %s\n' % ( command.rstrip())

  log_text += 'Version commands:\n'
  for (command, result) in zip(args.version_command_list, out_text):
    log_text += '  Cmd: %s\n' % (textwrap.indent(command.rstrip(), ' ' * 7).lstrip())
    log_text += '  Out: %s\n' % (textwrap.indent(result.rstrip(), ' ' * 7).lstrip())

  if(args.file_name != None):
    log_text += 'Additional logged information:\n'
    ifp = open(args.file_name, 'r')
    log_text += '  %s\n' % (textwrap.indent(ifp.read().rstrip(), ' ' * 2).lstrip())

  log_text += '\n'

  # Write to log file.
  if(args.output_file == None):
    log_file_name = args.sample_name + '.' + args.start_time.replace(':', '') + '.' + args.process_name + '.' + args.nextflow_run_name + '.' + work_id + '.log'
  else:
    log_file_name = args.output_file
  if(args.output_directory != None):
    log_file_name = os.path.join(args.output_directory, log_file_name)

  ofp = open(log_file_name, 'w')
  ofp.write(log_text)
  ofp.close()

