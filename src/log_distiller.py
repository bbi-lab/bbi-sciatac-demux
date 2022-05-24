#!/usr/bin/env python3

#
# Read log files written by pipeline_logger.py and distill
# the information into a single file with the requested
# format.
#
# The official repository for the program is
#
#   https://github.com/bbi-lab/bbi-sciatac-demux/src/log_distiller.py
#
# The program is distributed from there to
#
#   https://github.com/bbi-lab/bbi-sciatac-analyze/src/log_distiller.py
#

import os
import sys
import subprocess
import argparse
import textwrap
import json

LOG_DISTILLER_VERSION = '20220502a'


# Sample ordering lists. In fact these are non-sample-based processes such as
# lane processing and processes used by all samples.
sample_order_demux = ['lane', 'dashboard', 'all_samples']
sample_order_analyze = ['pipeline', 'genome', 'peaks', 'plots', 'dashboard']


# Process ordering lists are used to put process log files in processing order.
# The process names are in the order (or as close as possible) in which they are
# expected to be executed. Optional processes should be included as well.
process_order_demux = ['make_fake_sample_sheet', 'bcl2fastq', 'barcode_correct',
                       'adapter_trimming', 'fastqc_lanes', 'fastqc_samples', 'demux_dash']
process_order_analyze = ['sortTssBedProcess', 'sortChromosomeSizeProcess', 'sortPeakFileProcess',
                         'runAlignProcess', 'mergeBamsProcess', 'mergeMitoBamsProcess',
                         'getUniqueFragmentsProcess', 'getUniqueFragmentsMitoProcess',
                         'combineReadCountsProcess', 'callPeaksProcess',
                         'mergePeaksByGroupProcess', 'mergePeaksByFileProcess',
                         'makeWindowedGenomeIntervalsProcess', 'makePromoterSumIntervalsProcess',
                         'makeMergedPeakRegionCountsProcess', 'makeTssRegionCountsProcess',
                         'makeCountReportsProcess', 'callCellsProcess',
                         'getPerBaseCoverageTssProcess', 'makePeakMatrixProcess',
                         'makeWindowMatrixProcess', 'makePromoterMatrixProcess',
                         'summarizeCellCallsProcess', 'makeGenomeBrowserFilesProcess',
                         'getBandingScoresProcess', 'callMotifsProcess',
                         'makeMotifMatrixProcess', 'makeReducedDimensionMatrixProcess',
                         'experimentDashboardProcess', 'makeMergedPlotFilesProcess']


# Given a filename as a string, return a dictionary with
# the process name, sample name, run time, nextflow run id,
# and nextflow work directory hash.
# Split string on '.' then reassemble the sample name as
# required. An example filename is
#
#   ATAC.BAT1-RUN001_L001.20220420184701.adapter_trimming.crazy_lovelace.94_a7c1cc.json
#
# There is a minimum of 6 parts after the initial
# split on '.'
def filename_tokenizer(filename):
  parts  = filename.split('.')
  nparts = len(parts)
  if(nparts < 6):
    return({'is_logfile': False})

  is_logfile      = parts[-1] == 'json'
  work_directory  = parts[-2]
  pipeline_runid  = parts[-3]
  process_name    = parts[-4]
  process_runtime = int(parts[-5])
  sample_name_nparts = nparts - 5
  sample_name = '.'.join(parts[0:sample_name_nparts])
  filename_dict = {'is_logfile': is_logfile,
                   'sample_name': sample_name,
                   'pipeline_runid': pipeline_runid,
                   'process_name': process_name,
                   'process_runtime': process_runtime,
                   'work_directory': work_directory} 
  return(filename_dict)


# Make a list of process-level log filenames.
def make_log_file_list(args):
  process_list = args.processes
  # Store names of process-level log files in log_file_list.
  log_file_list = []
  if(args.log_file != None):
    log_directory, file_name = os.path.split(args.log_file)
    if(log_directory == ''):
      log_directory = '.'
    token_list = filename_tokenizer(file_name)
    if(not token_list['is_logfile']):
      print('Error: file \'%s\' is not a valid log file.' % (file_name))
      sys.exit(-1)
    log_file_list.append(file_name)
  elif(args.log_directory != None):
    log_directory = args.log_directory
    file_list = os.listdir(log_directory)
    # Check that the name splits into the expected number of
    # parts and has the .json suffix.
    for file_name in file_list:
      file_path = os.path.join(log_directory, file_name)
      # Ignore directories.
      if(not os.path.isfile(file_path)):
        continue
      token_list = filename_tokenizer(file_name)
      # Ignore files that are not process-level logs.
      if(not token_list['is_logfile']):
        continue
      if(not process_list is None and not token_list['process_name'] in process_list):
        continue
      if(args.samples is None or token_list['sample_name'] in args.samples): 
        log_file_list.append(file_name)
  return(log_directory, log_file_list)


# Make a dictionary with [sample][process][runid]
# where the runid has the most recent time stamp.
# These are not in the desired order; that is, the
# order in which the processing occurred, or close
# to it. The sample and process names and runids
# are taken from the log file names.
def make_sample_process_runid_dict(log_file_list):
  sample_process_runid_dict = {}
  for file_name in log_file_list:
    token_list = filename_tokenizer(file_name)
    # Store sample_name, process_name, and run_time
    sample_name     = token_list['sample_name']
    process_name    = token_list['process_name']
    pipeline_runid  = token_list['pipeline_runid']
    process_runtime = token_list['process_runtime']
    # Build sample_process_runid_dict.
    process_dict = sample_process_runid_dict.setdefault(sample_name, {})
    process_list = process_dict.setdefault(process_name, [None, None])
    if(process_list[1] is None or process_runtime > process_list[1]):
      sample_process_runid_dict[sample_name][process_name][0] = pipeline_runid
      sample_process_runid_dict[sample_name][process_name][1] = process_runtime
  return(sample_process_runid_dict)


# Make an ordered list of sample names used in this run.
def make_sample_ordered_list(sample_process_runid_dict, sample_order):
  # Make a list of ordered samples.
  sample_ordered_list = []
  # sample_flag_dict keys are sample names extracted from file names.
  # We pull out first the 'sample names' common to all runs, which
  # are listed in order in sample_order, and then append the run-
  # specific sample names.
  sample_flag_dict = dict.fromkeys(sample_process_runid_dict, 0)
  for a_sample in sample_order:
    if(a_sample in sample_flag_dict and sample_flag_dict[a_sample] == 0):
      sample_ordered_list.append(a_sample)
      sample_flag_dict[a_sample] = 1
  sample_ordered_list_tmp = []
  for a_sample in sample_flag_dict.keys():
    if(sample_flag_dict[a_sample] == 0):
      sample_ordered_list_tmp.append(a_sample)
  sample_ordered_list_tmp.sort()
  sample_ordered_list = sample_ordered_list + sample_ordered_list_tmp
  return(sample_ordered_list)


# Make a dictionary of samples with ordered process lists as values.
def make_sample_process_ordered_dict(sample_ordered_list, sample_process_runid_dict, process_order):
  # Make a dictionary keyed by sample whose values are ordered
  # lists of processes.
  sample_process_ordered_dict = {}
  for a_sample in sample_ordered_list:
    a_process_dict = sample_process_runid_dict[a_sample]
    process_ordered_list = []
    process_flag_dict = dict.fromkeys(a_process_dict, 0)
    for a_process in process_order:
      if(a_process in process_flag_dict and process_flag_dict[a_process] == 0):
        process_ordered_list.append(a_process)
        process_flag_dict[a_process] = 1
    process_ordered_list_tmp = []
    for a_process in process_flag_dict:
      if(process_flag_dict[a_process] == 0):
        process_ordered_list_tmp.append(a_process)
    process_ordered_list = process_ordered_list + process_ordered_list_tmp
    sample_process_ordered_dict[a_sample] = process_ordered_list
  return(sample_process_ordered_dict)


# Select process-level files and sort them.
def select_log_files(args):
  if(args.processing_pipeline == 'atac_demux'):
    sample_order = sample_order_demux
    process_order = process_order_demux
  elif(args.processing_pipeline == 'atac_analyze'):
    sample_order = sample_order_analyze
    process_order = process_order_analyze

  log_directory, log_file_list = make_log_file_list(args)

  sample_process_runid_dict    = make_sample_process_runid_dict(log_file_list)
  sample_ordered_list          = make_sample_ordered_list(sample_process_runid_dict, sample_order)
  sample_process_ordered_dict  = make_sample_process_ordered_dict(sample_ordered_list, sample_process_runid_dict, process_order)

  # Set up a dictionary by sample of dictionaries by process of lists
  # of corresponding log files.
  sample_process_runid_files = dict.fromkeys(sample_ordered_list)
  for a_sample in sample_process_runid_files.keys():
    sample_process_runid_files[a_sample] = dict.fromkeys(sample_process_ordered_dict[a_sample])
    for a_process in sample_process_runid_files[a_sample].keys():
      sample_process_runid_files[a_sample][a_process] = []
  
  # Populate the sample_process_runid_files list.
  for a_file in log_file_list:
    token_list = filename_tokenizer(a_file)
    sample_name    = token_list['sample_name']
    process_name   = token_list['process_name']
    pipeline_runid = token_list['pipeline_runid']

    if(pipeline_runid == sample_process_runid_dict[sample_name][process_name][0]):
      a_path = os.path.join(log_directory, a_file)
      sample_process_runid_files[sample_name][process_name].append([a_file, a_path])

  return({'sample_ordered_list': sample_ordered_list, 'sample_process_ordered_dict': sample_process_ordered_dict, 'sample_process_runid_files': sample_process_runid_files})


# Read a process-level log file and return the contents
# as a dictionary. The file contents are in JSON format
# so this function simply reads the file into a dictionary.
def read_logfile(file_name):
  file_contents_dict = json.load(open(file_name))
  return(file_contents_dict)


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


# Indent text block but not the leading line.
def indent_text_block(text_block, spacer):
  parts = text_block.split('\n')
  indent = '\n' + spacer
  out_text_block = indent.join(parts)
  return(out_text_block)
  

def process_log_files(log_file_dict, output_file=sys.stdout):
  for a_sample in log_file_dict['sample_ordered_list']:
    print('==== sample: %s ====\n' % (a_sample), file=output_file)
    for a_process in log_file_dict['sample_process_ordered_dict'][a_sample]:
      print('  ** process: %s **\n' % (a_process), file=output_file)
      for a_list in log_file_dict['sample_process_runid_files'][a_sample][a_process]:
        file_contents_dict = read_logfile(a_list[1])
        indent_level = 1
        write_data_dict(file_contents_dict, indent_level, output_file)
      print(file=output_file)
    print(file=output_file)


if __name__ == '__main__':
  parser = argparse.ArgumentParser('Script to distill process-level log files made by pipeline_logger.py into a single report.')
  parser.add_argument('-p', '--processing_pipeline', required=True, help='The pipeline that made the log files. Either "atac_demux" or "atac_analyze".')
  parser.add_argument('-t', '--log_file', required=False, default=None, help='The name of a process-level log file (used for development).')
  parser.add_argument('-d', '--log_directory', required=False, default=None, help='The name of the directory that has the process-level log files.')
  parser.add_argument('-l', '--log_filter_level', required=False, default='all', help='The amount of information to include in the output file. Default is all.')
  parser.add_argument('-s', '--samples', required=False, nargs='*', default=None, help='The samples to report in the output file. Default is all.')
  parser.add_argument('-r', '--processes', required=False, nargs='*', default=None, help='The processes to report in the output file. Default is all.')
  parser.add_argument('-o', '--output_file', required=False, default=None, help='Optional output filename. Default to stdout.')
  args = parser.parse_args()

if(not args.log_file is None and not args.log_directory is None):
  print('Error: you must specify either -f or -d but not both.')
  sys.exit(-1)
elif(args.log_file is None and args.log_directory is None):
  print('Error: you must specify either -f or -d.')
  sys.exit(-1)

if(args.output_file is None):
  output_file = sys.stdout
else:
  output_file = open(args.output_file, 'w+')

log_file_dict = select_log_files(args)
process_log_files(log_file_dict, output_file)

if(not args.output_file is None):
  output_file.close()


