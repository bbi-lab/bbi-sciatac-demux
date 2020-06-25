from __future__ import print_function
from __future__ import division
import argparse
import subprocess
import sys
import os
import gzip
import io
import itertools
import time
import json
import re
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

def indexsplitter(indexrange):
	if len(indexrange) < 3:
		indexout = [int(indexrange)-1]
	elif "-" in indexrange or ',' in indexrange:
		range_list = [x for x in indexrange.split(",")]
		indexout = []
		for myrange in range_list:
			index_range = myrange.split('-')
			
			if len(index_range) == 1:
				start = int(index_range[0]) - 1
				end = start + 1
			elif len(index_range) == 2:
				start = int(index_range[0]) - 1
				end = int(index_range[1])
			else:
				raise ValueError('Invalid index range %s' % myrange)

			indexout.extend(range(start, end))
	else:
		raise ValueError('Invalid format for index range: %s' % indexrange)
	return indexout



def required_index(a):
    """
    Helper function to take a list of index lists and return whether it needs to be included as an index in demultiplexing.
    """
    return len(set(tuple(a_i) for a_i in a)) != 1
    
def get_sample_lookup(samplesheet, pcri7, tagi7, tagi5, pcri5):
    """
    Takes a samplesheet file handle and returns a list of True/False indicating whether index was used in lookup and a lookup of tuples of index combinations to sample names.

    Args:
        samplesheet (file handle): opened file handle to samplesheet
        tagi5 (list): list of tag/lig i5 indices
        pcri5 (list): list of pcr i5 indices
        pcri7 (list): list of pcr i7 indices
        tagi7 (list): list of tag/lig i7 indices

    Returns:
        (list of bool, dict of tuple to sample name, dict of valid tagmentation index pair tuples, dict of valid PCR index pair tuples): [use_pcri7, use_tagi7, use_tagi5, use_pcri5] where each entry indicates usage of that barcode in indexing and lookup of (tagi7,tagi5), for example.

    """
    pcri7_indices = []
    tagi7_indices = []
    tagi5_indices = []
    pcri5_indices = []
    samples = []
    for line in samplesheet:
        if line.startswith('sample_id\tranges'):
            continue

        entries = line.strip().split()
        sample, indices, genome = entries
        indexsplit = indices.split(':')

        samples.append(sample)
        pcri7_indices.append(indexsplitter(indexsplit[1]))
        tagi7_indices.append(indexsplitter(indexsplit[0]))
        tagi5_indices.append(indexsplitter(indexsplit[3]))
        pcri5_indices.append(indexsplitter(indexsplit[2]))

    use_pcri7 = required_index(pcri7_indices)
    use_tagi7 = required_index(tagi7_indices)
    use_tagi5 = required_index(tagi5_indices)
    use_pcri5 = required_index(pcri5_indices)

    index_mask = [use_pcri7, use_tagi7, use_tagi5, use_pcri5]
    index_whitelists = [pcri7, tagi7, tagi5, pcri5]
    index_lists = [pcri7_indices, tagi7_indices, tagi5_indices, pcri5_indices]
    
    sample_lookup_table = {}
    for sample_index,sample in enumerate(samples):
        indices_to_use = []
        for i,(use_index,index_list) in enumerate(zip(index_mask, index_lists)):
            if use_index:
                indices_to_use.append([index_whitelists[i][index_i] for index_i in index_list[sample_index]])

        for combination in list(itertools.product(*indices_to_use)):
            sample_lookup_table[combination] = sample
    
    pcr_pairs = {}
    for (pcri7_list, pcri5_list) in zip( pcri7_indices, pcri5_indices):
      for combination in itertools.product(pcri7_list, pcri5_list):
        pcr_pairs.setdefault( (combination[0]-1, combination[1]-1), 1 )

    tag_pairs = {}
    for (tagi7_list, tagi5_list) in zip( tagi7_indices, tagi5_indices):
      for combination in itertools.product(tagi7_list, tagi5_list):
        tag_pairs.setdefault( (combination[0]-1, combination[1]-1), 1 )

    return index_mask, sample_lookup_table, tag_pairs, pcr_pairs


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='A program to fix erroneous barcodes in scATAC data.')
    parser.add_argument('-1', '--input1', nargs='?', type=argparse.FileType('r'), default=sys.stdin, required=True, help='Text piped in from stdin for R1.')
    parser.add_argument('-2', '--input2', nargs='?', type=argparse.FileType('r'), default=sys.stdin, required=True, help='Text piped in from stdin for R2.')
    parser.add_argument('--filename', required=True, help='The R1 file name.')
    parser.add_argument('--samplesheet', required=True, help='Samplesheet describing the layout of the samples.')
    parser.add_argument('--out_dir', required=True, help='Output directory.')
    parser.add_argument('--stats_out', required=True, help='write JSON file with output stats about processed reads and correction rates (0=no;1=yes.')
    parser.add_argument('--two_level_indexed_tn5', action='store_true', help='Flag to run assuming that the library is a two-level indexed-TN5 sample.')
    parser.add_argument('--wells_384', action='store_true', help='Flag to run assuming That the known barcode set is the 384 well set.')
    parser.add_argument('--well_ids', action='store_true', help='Flag to output cell IDs that are composed of well IDs rather than the actual sequences.')
    parser.add_argument('-X', '--nextseq', help='NextSeq run indicator', dest='nextseq', action="store_true")
    args = parser.parse_args()

    # Note: args.filename may include (leading) directory information.
    lane_num = args.filename
    lane_num = re.sub('.*Undetermined_S0_L','L',lane_num)
    lane_num = lane_num.replace('_R1_001.fastq.gz', '')

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

    # Build up sample mapping from indices to samples
    index_mask, sample_lookup, tag_pairs, pcr_pairs  = get_sample_lookup(open(args.samplesheet), pcri7, tagi7, tagi5, pcri5)

    tagmentation_i7_whitelist = tagi7
    pcr_i7_whitelist = pcri7

    if not args.two_level_indexed_tn5:

        if args.nextseq:
            p5_pcr_rc_map = {reverse_complement(k):k for k in pcri5}
            p5_tagmentation_rc_map = {reverse_complement(k):k for k in tagi5}

            pcr_i5_whitelist = set([reverse_complement(x) for x in pcri5])
            tagmentation_i5_whitelist = set([reverse_complement(x) for x in tagi5])
        else:
            pcr_i5_whitelist = pcri5
            tagmentation_i5_whitelist = tagi5

    else:
        tagmentation_i7_whitelist = bc.nex_i7_two_level_indexed_tn5
        pcr_i7_whitelist = bc.pcr_i7_two_level_indexed_tn5

        if args.nextseq:
            p5_pcr_rc_map = {reverse_complement(k):k for k in bc.pcr_i5_two_level_indexed_tn5}
            p5_tagmentation_rc_map = {reverse_complement(k):k for k in bc.nex_i5_two_level_indexed_tn5}

            pcr_i5_whitelist = set([reverse_complement(x) for x in bc.pcr_i5_two_level_indexed_tn5])
            tagmentation_i5_whitelist = set([reverse_complement(x) for x in bc.nex_i5_two_level_indexed_tn5])
        else:
            pcr_i5_whitelist = bc.pcr_i5_two_level_indexed_tn5
            tagmentation_i5_whitelist = bc.nex_i5_two_level_indexed_tn5

    tagi7_to_index = barcode_to_well.barcode_index_dict( tagmentation_i7_whitelist )
    pcri7_to_index = barcode_to_well.barcode_index_dict( pcr_i7_whitelist )
    pcri5_to_index = barcode_to_well.barcode_index_dict( pcr_i5_whitelist )
    tagi5_to_index = barcode_to_well.barcode_index_dict( tagmentation_i5_whitelist )

    tagmentation_i7_correction_map = construct_mismatch_to_whitelist_map(tagmentation_i7_whitelist, 2)
    pcr_i7_correction_map = construct_mismatch_to_whitelist_map(pcr_i7_whitelist, 2)
    pcr_i5_correction_map = construct_mismatch_to_whitelist_map(pcr_i5_whitelist, 2)
    tagmentation_i5_correction_map = construct_mismatch_to_whitelist_map(tagmentation_i5_whitelist, 2)

    # Set up all input/output files
    if( not os.path.exists(args.out_dir)):
        os.mkdir(args.out_dir)
    
    output_files = {}
    for sample in list(set(sample_lookup.values())):
        output_file_1 = os.path.join(args.out_dir, '%s-RUN001_%s_R1.fastq' % (sample, lane_num))
        output_file_2 = os.path.join(args.out_dir, '%s-RUN001_%s_R2.fastq' % (sample, lane_num))
        output_files[sample] = {}
        output_files[sample]['r1'] = open(output_file_1, 'w')
        output_files[sample]['r1_name'] = output_file_1
        output_files[sample]['r2'] = open(output_file_2, 'w')
        output_files[sample]['r2_name'] = output_file_2

    
    output_file_stats_json = os.path.join(args.out_dir, 'RUN001_%s.stats.json' % (lane_num))
    output_file_counts_barcodes_csv = os.path.join(args.out_dir, 'RUN001_%s.barcode_counts.csv' % (lane_num))
    output_file_counts_indexes_csv = os.path.join(args.out_dir, 'RUN001_%s.index_counts.csv' % (lane_num))

    if1 = FastqGeneralIterator(args.input1)
    if2 = FastqGeneralIterator(args.input2)

    totreads = 0
    total_not_specified_in_samplesheet = 0
    validreads = {}
    validreads['pcr_i5'] = 0
    validreads['pcr_i7'] = 0
    validreads['pcr'] = 0
    validreads['pcr_match'] = 0
    validreads['tagmentation_i5'] = 0
    validreads['tagmentation_i7'] = 0
    validreads['tagmentation'] = 0
    validreads['tagmentation_match'] = 0
    validreads['all_barcodes'] = 0

    start = time.time()

    barcodes_string_count = {}
    tagmentation_i7_count = [0] * 384
    pcr_i7_count = [0] * 384
    pcr_i5_count = [0] * 384
    tagmentation_i5_count = [0] * 384
 
    for (r1_name, r1_seq, r1_qual),(r2_name, r2_seq, r2_qual) in zip(if1, if2):

        totreads += 1

        # Get barcodes and correct
        tagmentation_i7_seq, pcr_i7_seq, pcr_i5_seq, tagmentation_i5_seq = get_barcode_seqs(r1_name, args.nextseq, args.two_level_indexed_tn5)
        tagi7_in = tagmentation_i7_seq
        pcri7_in = pcr_i7_seq
        pcri5_in = pcr_i5_seq
        tagi5_in = tagmentation_i5_seq
        tagmentation_i7_seq = correct_barcode(tagmentation_i7_seq, tagmentation_i7_correction_map)
        pcr_i7_seq = correct_barcode(pcr_i7_seq, pcr_i7_correction_map)
        pcr_i5_seq = correct_barcode(pcr_i5_seq, pcr_i5_correction_map)
        tagmentation_i5_seq = correct_barcode(tagmentation_i5_seq, tagmentation_i5_correction_map)
      
        # Skip invalid reads and track valid read count for error checking
        if tagmentation_i7_seq is not None:
            validreads['tagmentation_i7'] += 1
            tagmentation_i7_index = tagi7_to_index.get(tagmentation_i7_seq)
            if tagmentation_i7_index is None:
                print( 'tagmentation i7 index is none for sequence %s (%s)' % (tagmentation_i7_seq,tagi7_in) )
            tagmentation_i7_count[tagmentation_i7_index] += 1
        if tagmentation_i5_seq is not None:
            validreads['tagmentation_i5'] += 1
            tagmentation_i5_index = tagi5_to_index.get(tagmentation_i5_seq)
            if tagmentation_i5_index is None:
                print( 'tagmentation i5 index is none for sequence % (%s)s' % (tagmentation_i5_seq,tagi5_in) )
            tagmentation_i5_count[tagmentation_i5_index] += 1
        if tagmentation_i7_seq is not None and tagmentation_i5_seq is not None:
            validreads['tagmentation'] += 1
            if tag_pairs[(tagmentation_i7_index, tagmentation_i5_index)] is not None:
                validreads['tagmentation_match'] += 1
        if pcr_i7_seq is not None:
            validreads['pcr_i7'] += 1
            pcr_i7_index = pcri7_to_index.get(pcr_i7_seq)
            if pcr_i7_index is None:
                print( 'pcr i7 index is none for sequence %s (%s)' % (pcr_i7_seq,pcri7_in) )
            pcr_i7_count[pcr_i7_index] += 1
        if pcr_i5_seq is not None:
            validreads['pcr_i5'] += 1
            pcr_i5_index = pcri5_to_index.get(pcr_i5_seq)
            if pcr_i5_index is None:
                print( 'pcr i5 index is none for sequence %s (%s)' % (pcr_i5_seq,pcri5_in) )
            pcr_i5_count[pcr_i5_index] += 1
        if pcr_i7_seq is not None and pcr_i5_seq is not None:
            validreads['pcr'] += 1
            if pcr_pairs[(pcr_i7_index, pcr_i5_index)] is not None:
                validreads['pcr_match'] += 1
        
        if tagmentation_i7_seq is None or pcr_i7_seq is None or pcr_i5_seq is None or tagmentation_i5_seq is None:
            continue

        validreads['all_barcodes'] += 1

        # Map back to original whitelist if on nextseq so barcodes are always same on every sequencer
        if args.nextseq:
            pcr_i5_seq = p5_pcr_rc_map[pcr_i5_seq]
            tagmentation_i5_seq = p5_tagmentation_rc_map[tagmentation_i5_seq]

        # Convert to well IDs if requested
        # Note that for two level Tn5 barcode well comes first then PCR,
        # for three-level will be Tn5 N7, Tn5 N5, PCR WELL ID
        if args.two_level_indexed_tn5:
            barcodes_string = barcode_to_well.get_two_level_barcode_string(tagmentation_i7_seq, pcr_i7_seq, pcr_i5_seq, tagmentation_i5_seq, bc.nex_two_level_indexed_tn5_to_well, bc.pcr_two_level_indexed_tn5_to_well, args.well_ids)
        else:
            barcodes_string = barcode_to_well.get_barcode_string(tagmentation_i7_seq, pcr_i7_seq, pcr_i5_seq, tagmentation_i5_seq, lig_i7_to_well, lig_i5_to_well, pcr_to_well, args.well_ids)

        sample_index = tuple([index for use_index,index in zip(index_mask, [pcr_i7_seq, tagmentation_i7_seq, tagmentation_i5_seq, pcr_i5_seq]) if use_index])
        sample = sample_lookup.get(sample_index, None)

        if not sample:
            total_not_specified_in_samplesheet += 1
        else:
            output_files[sample]['r1'].write(''.join(['@', barcodes_string, ':', str(totreads), '\n', r1_seq, '\n+\n', r1_qual, '\n']))
            output_files[sample]['r2'].write(''.join(['@', barcodes_string, ':', str(totreads), '\n', r2_seq, '\n+\n', r2_qual, '\n']))

        barcodes_string_count[barcodes_string] = barcodes_string_count.setdefault(barcodes_string, 0) + 1

    if totreads == 0:
        raise ValueError('No reads found in fastq input.')

    # Output basic stats
    # Note: write counts; not fractions
#    for stat in validreads:
#        validreads[stat] = validreads[stat] / float(totreads)
        
    validreads['total_input_reads'] = totreads
    validreads['total_not_specified_in_samplesheet'] = total_not_specified_in_samplesheet
#    validreads['fraction_not_specified_in_samplesheet'] = total_not_specified_in_samplesheet / totreads

    with open(output_file_stats_json, 'wt') as f:
        f.write(json.dumps(validreads, indent=4))

    with open(output_file_counts_barcodes_csv,'wt') as f:
        for barcodes_string in barcodes_string_count.keys():
            btmp = barcodes_string if (type( barcodes_string ) == str) else 'None'
            f.write(''.join([lane_num, ',', btmp, ',', str(barcodes_string_count[barcodes_string]),'\n']))
            
    with open(output_file_counts_indexes_csv,'wt') as f:
        for i, counts in enumerate(zip(tagmentation_i7_count, pcr_i7_count, pcr_i5_count, tagmentation_i5_count), start = 0):
            f.write(''.join([i,',',get_well_id_384_to_96(i, row_ordered = TRUE),',',counts[0],',',counts[1],',',get_well_id_384_to_96(i, row_ordered = FALSE),',',counts[2],',',counts[3],'\n']))
            
    # Error checking and compress output
    if validreads['all_barcodes'] < 0.05:
        raise ValueError('Warning, you had less than 5 percent of all reads pass index correction. Something may have gone wrong here w.r.t. index sets or the expected library configuration not matching the data...')

    print('Done correcting barcodes in %s minutes. Starting compression...' % ((time.time() - start) / 60.0))
    start = time.time()
    for sample in output_files:
        output_files[sample]['r1'].close()
        output_files[sample]['r2'].close()
        subprocess.check_call('pigz %s' % output_files[sample]['r1_name'], shell=True)
        subprocess.check_call('pigz %s' % output_files[sample]['r2_name'], shell=True)
    print('Done compressing with pigz in %s minutes.' % ((time.time() - start) / 60.0))
