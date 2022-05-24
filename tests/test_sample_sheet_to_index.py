#!/usr/bin/env python3


# Program: test_sample_sheet_to_index.py
# Purpose: explore and test the get_sample_lookup() function, which is in the file
#          bbi-sciatac-demux/src/barcode_correct_sciatac.py. This function creates
#          the dictionary of oligo sequence tuples that is used for demuxing reads.
#


# Notes:
#   o  if there is less than two samples the use_mask for all indices is False
#   o  if the indices are the same across all samples for an index type; e.g. p7 or p5,
#      then that index set is dropped; that is, it is not in the sample_lookup sequences
#      for all of the samples. Often, the n7, p7, and p5 are the same for all samples
#      and so the sample lookup list consists only of the and n5 sequences.
#   o  this program calles get_sample_lookup() with the no_mask parameter set to True
#      in order to simplify formatting the output table.


import argparse
import sys
import os
import gzip
import io
import itertools
import json
import re
import collections
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import barcode_to_well
import barcode_constants as bc
import barcode_correct_sciatac as bcc


#
# Notes:
#   o  test_get_sample_lookup tests the function get_sample_lookup() in
#      the file barcode_correct_sciatac.py
#   o  copy the following files from bbi-sciatac-demux/src into
#      this directory
#         o  barcode_constants.py
#         o  barcode_correct_sciatac.py
#         o  barcode_to_well.py
#
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='A program to check sample sheet indices.')
    parser.add_argument('--samplesheet', required=True, help='Samplesheet describing the layout of the samples.')
    parser.add_argument('--two_level_indexed_tn5', action='store_true', help='Flag to run assuming that the library is a two-level indexed-TN5 sample.')
    parser.add_argument('--wells_384', action='store_true', help='Flag to run assuming that the known barcode set is the 384 well set.')
    parser.add_argument('-X', '--nextseq', help='NextSeq run indicator', dest='nextseq', action="store_true")
    parser.add_argument('--well_ids', action='store_true', help='Flag to output cell IDs that are composed of well IDs rather than the actual sequences.')
    parser.add_argument('--no_mask', action='store_true', help='Use all four barcodes. By default, do not use (mask out) a barcode(s) when all samples have the same index set (flag).')
 
    args = parser.parse_args()

    if args.two_level_indexed_tn5 and args.wells_384:
        raise ValueError('There is no 384 well barcode set for indexed Tn5, may not specify both --two_level_indexed_tn5 and --wells_384.')

    # Set up the right index set depending on the indices
    if args.two_level_indexed_tn5:
        tagi7 = bc.nex_i7_two_level_indexed_tn5_list
        pcri7 = bc.pcr_i7_two_level_indexed_tn5_list
        pcri5 = bc.pcr_i5_two_level_indexed_tn5_list
        tagi5 = bc.nex_i5_two_level_indexed_tn5_list
    else:
        if args.wells_384:
            tagi7 = bc.lig_i7_list_384
            pcri7 = bc.pcr_i7_list_384
            pcri5 = bc.pcr_i5_list_384
            tagi5 = bc.lig_i5_list_384
            lig_i7_to_well = bc.lig_i7_to_well_384
            lig_i5_to_well = bc.lig_i5_to_well_384
            pcr_to_well = bc.pcr_to_well_384
        else:
            tagi7 = bc.lig_i7_list
            pcri7 = bc.pcr_i7_list
            pcri5 = bc.pcr_i5_list
            tagi5 = bc.lig_i5_list
            lig_i7_to_well = bc.lig_i7_to_well
            lig_i5_to_well = bc.lig_i5_to_well
            pcr_to_well = bc.pcr_to_well

    tagi7_index_dict = barcode_to_well.barcode_index_dict(tagi7)
    pcri7_index_dict = barcode_to_well.barcode_index_dict(pcri7)
    pcri5_index_dict = barcode_to_well.barcode_index_dict(pcri5)
    tagi5_index_dict = barcode_to_well.barcode_index_dict(tagi5)

    # Build up sample mapping from indices to samples
    index_mask, sample_lookup, tagi5_sample_list, tag_pairs_counts, pcr_pairs_counts, index_flags = bcc.get_sample_lookup(open(args.samplesheet), True, pcri7, tagi7, tagi5, pcri5)

    print('sample  i7_index  p7_index  p5_index  i5_index')
    for sample_tuple, sample in sample_lookup.items():
        print('%s    %d %d %d %d' % ( sample, tagi7_index_dict[sample_tuple[1]]+1, pcri7_index_dict[sample_tuple[0]]+1, pcri5_index_dict[sample_tuple[3]]+1, tagi5_index_dict[sample_tuple[2]]+1 ) )



