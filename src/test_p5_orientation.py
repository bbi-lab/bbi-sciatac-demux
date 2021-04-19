#!/usr/bin/env python3

#
# Program: test_p5_orientation.py
# Application: Ilumina sci-ATAC-seq sequencing.
# Purpose: try to determine the P5 sequence
#          orientation by reading a small set of
#          reads from a bcl2fastq output fastq file.
# Notes:
#   o  this code is from Andrew Hill's
#      barcode_correct_sciatac.py program.
#   o  my hope is to use this program to
#      avoid the need to infer the P5 sequence
#      orientation from the RunParameters.xml
#      file. The immediate motivation is that
#      the P5 orientation for the NovaSeq
#      v1.5 reagents appears to depend on
#      whether the experiment is RNA-seq or
#      ATAC-seq, in addition to the information
#      in RunParameters.xml. I think that using
#      the correct P5 orientation is important
#      for NovaSeq runs in order to avoid wasting
#      time producing empty fastq files. (Although
#      it appears that making empty demuxed
#      fastq files takes only 3 or 4 hours.)
#   o  I am concerned about possible errors
#      in the orientation determination due to
#      things like bad library preparations.
#      I would like to test this program on a
#      variety of data sets in order to get a
#      sense of the range of values that I can
#      expect, and consider to be valid.
#

import argparse
import io
import itertools
import json
import re
import gzip
import collections
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import barcode_to_well
import barcode_constants as bc


def correct_barcode(barcode, mismatch_map):
    """
    Correct an observed raw barcode to one of a list of whitelists of mismatches.
    Args:
            barcode (string): barcode sequence to be corrected
            mismatch_map (list of dict dict): list of dict of mismatched sequences to real sequences
    Returns:
            string: corrected barcodes or None if barcode not correctable.
    """
    for mismatch_whitelist in mismatch_map:
        corrected = mismatch_whitelist.get(barcode, None)

        if corrected:
            return corrected

    return None


def generate_mismatches(sequence, num_mismatches, allow_n=True):
    """
    Generate a list of mismatched sequences to a given sequence. Must only contain ATGC.
    This is heavily based on a biostars answer.
    Args:
        sequence (str): The sequence must contain only A, T, G, and C
        num_mismatches (int): number of mismatches to generate sequences for
        allow_n (bool): True to allow N bases and False if not
    Yield:
    """
    letters = 'ACGT'

    if allow_n:
        letters += 'N'

    sequence = sequence.upper()
    mismatches = []

    # generate combinations of num_mismatch indices into the sequence
    for locs in itertools.combinations(range(len(sequence)), num_mismatches):
        # convert sequence to list of characters
        sequence_list = [[char] for char in sequence]
        for loc in locs:
            orig_char = sequence[loc]
            # replace target character with list of mismatches
            sequence_list[loc] = [l for l in letters if l != orig_char]

        # expand lists as cartesian product (and convert list of characters to string)
        for poss in itertools.product(*sequence_list):
            mismatches.append(''.join(poss))

    return mismatches


def construct_mismatch_to_whitelist_map(whitelist, edit_distance, allow_n=True):
    """
    Constructs a precomputed set of all mismatches within a specified edit distance and the barcode whitelist.
    Args:
        whitelist (set of str): set of whitelist sequences
        edit_distance (int): max edit distance to consider
        allow_n (bool): True to allow N bases and False if not
    Returns:
        dict: mapping of mismatched sequences to their whitelist sequences
    """

    mismatch_to_whitelist_map = [None] * (edit_distance + 1)

    mismatch_to_whitelist_map[0] = {k: k for k in whitelist}

    conflicting_mismatches = []  # tracks conflicts where mismatches map to different sequences

    # Doesn't really matter as correction function will never see it,
    # but exclude any perfect matches to actual seqs by mismatches
    conflicting_mismatches.extend(list(whitelist))

    for mismatch_count in range(1, edit_distance + 1):
        mismatch_to_whitelist_map[mismatch_count] = {}

        for sequence in whitelist:
            sequence = sequence.upper()

            # Generate all possible mismatches in range
            mismatches = generate_mismatches(sequence, num_mismatches=mismatch_count, allow_n=allow_n)

            # Construct a mapping to the intended sequences
            for mismatch in mismatches:
                # Check for conflict with existing sequence and track if so
                if mismatch in mismatch_to_whitelist_map[mismatch_count]:
                    conflicting_mismatches.append(mismatch)
                mismatch_to_whitelist_map[mismatch_count][mismatch] = sequence

        # Go back and remove any conflicting mismatches
        for mismatch in set(conflicting_mismatches):
            if mismatch in mismatch_to_whitelist_map[mismatch_count]:
                del mismatch_to_whitelist_map[mismatch_count][mismatch]

    return mismatch_to_whitelist_map


def reverse_complement(x):
    complements = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    xrev = x[::-1]
    xrevcomp = ''.join([complements[z] for z in xrev])
    return xrevcomp


def get_barcode_seqs(r1_name, nextseq, two_level_indexed_tn5):
    """
    Extract the correct sequences from the R1 name.
    """
    # In 3LV runs, this is the the P5 + P7 index seq with a + in between
    # which is 20 + 1 + 20 (so have to skip "+" position)
    # Similar for two-level indexed Tn5, but Tn5 barcodes are 8bp
    if not two_level_indexed_tn5:
        barcodes = r1_name[-41:]

        tagmentation_i7_seq = barcodes[0:10]
        pcr_i7_seq = barcodes[10:20]

        if nextseq:
            pcr_i5_seq = barcodes[31:41]
            tagmentation_i5_seq = barcodes[21:31]
        else:
            pcr_i5_seq = barcodes[21:31]
            tagmentation_i5_seq = barcodes[31:41]
    else:
        barcodes = r1_name[-37:]

        tagmentation_i7_seq = barcodes[0:8]
        pcr_i7_seq = barcodes[8:18]

        if nextseq:
            pcr_i5_seq = barcodes[27:37]
            tagmentation_i5_seq = barcodes[19:27]
        else:
            pcr_i5_seq = barcodes[19:29]
            tagmentation_i5_seq = barcodes[29:37]

    return tagmentation_i7_seq, pcr_i7_seq, pcr_i5_seq, tagmentation_i5_seq


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='A program to fix erroneous barcodes in scATAC data.')
    parser.add_argument('--filename', required=True, help='The R1 file name.')
    parser.add_argument('--two_level_indexed_tn5', action='store_true', help='Flag to run assuming that the library is a two-level indexed-TN5 sample.')
    parser.add_argument('--wells_384', action='store_true', help='Flag to run assuming that the known barcode set is the 384 well set.')

    args = parser.parse_args()

    if args.two_level_indexed_tn5 and args.wells_384:
        raise ValueError('There is no 384 well barcode set for indexed Tn5, may not specify both --two_level_indexed_tn5 and --wells_384.')

    if args.two_level_indexed_tn5:
        pcri5 = bc.pcr_i5_two_level_indexed_tn5_list
    else:
        if args.wells_384:
            pcri5 = bc.pcr_i5_list_384
        else:
            pcri5 = bc.pcr_i5_list

    if not args.two_level_indexed_tn5:

            p5_pcr_rc_map = {reverse_complement(k):k for k in pcri5}
            pcr_i5_whitelist_rev = set([reverse_complement(x) for x in pcri5])
            pcr_i5_whitelist_fwd = pcri5
    else:
            p5_pcr_rc_map = {reverse_complement(k):k for k in bc.pcr_i5_two_level_indexed_tn5}
            pcr_i5_whitelist_rev = set([reverse_complement(x) for x in bc.pcr_i5_two_level_indexed_tn5])
            pcr_i5_whitelist_fwd = bc.pcr_i5_two_level_indexed_tn5

    pcr_i5_correction_map_fwd = construct_mismatch_to_whitelist_map(pcr_i5_whitelist_fwd, 2)
    pcr_i5_correction_map_rev = construct_mismatch_to_whitelist_map(pcr_i5_whitelist_rev, 2)

    fp = gzip.open(args.filename, 'rt')
    if1 = FastqGeneralIterator(fp)

    num_count = 10000
    beg_count = 1000000
    end_count = beg_count + num_count

    tot_reads = 0 
    num_pcr_i5_fwd = 0
    num_pcr_i5_rev = 0
    for (r1_name, r1_seq, r1_qual) in if1:

        tot_reads += 1

        if(tot_reads < beg_count):
            continue
        elif( tot_reads > end_count):
            break

        # Get barcodes and correct
        # Ignore unused index sequences.
        nextseq = False
        tagmentation_i7_seq_fwd, pcr_i7_seq_fwd, pcr_i5_seq_fwd, tagmentation_i5_seq_fwd = get_barcode_seqs(r1_name, nextseq, args.two_level_indexed_tn5)
        nextseq = True
        tagmentation_i7_seq_rev, pcr_i7_seq_rev, pcr_i5_seq_rev, tagmentation_i5_seq_rev = get_barcode_seqs(r1_name, nextseq, args.two_level_indexed_tn5)

        pcr_i5_seq_fwd = correct_barcode(pcr_i5_seq_fwd, pcr_i5_correction_map_fwd)
        pcr_i5_seq_rev = correct_barcode(pcr_i5_seq_rev, pcr_i5_correction_map_rev)
        if(pcr_i5_seq_fwd != None):
            num_pcr_i5_fwd += 1
        if(pcr_i5_seq_rev != None):
            num_pcr_i5_rev += 1

    fp.close()

    frac_fwd = float(num_pcr_i5_fwd) / float(num_count)
    frac_rev = float(num_pcr_i5_rev) / float(num_count)

    print( 'frac_fwd: %f' % ( frac_fwd ) )
    print( 'frac_rev: %f' % ( frac_rev ) )


