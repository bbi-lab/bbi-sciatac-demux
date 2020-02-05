from __future__ import print_function

#Script to generate a table of indexes and sample names
#For each sample name, input the nextera i7 index numbers,
#pcr i7 index numbers, pcr i5 index numbers and nextera i5
#index numbers in the format lig_i7:pcr_i7:pcr_i5:lig_i5.
#To assign multiple indexes, separate numbers between colons
#with commas (e.g. '1,3,5').  To input a range, separate
#start and end with '-' (e.g. '1-6').
import sys
import argparse
import barcode_to_well
import barcode_constants as bc

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


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Script to generate an index table for sci-ATAC runs.')
	parser.add_argument('indices', help='Index string. See overall pipeline docs for how these work.')
	parser.add_argument('name', help='A name for the sample who fit this index range.')
	parser.add_argument('--two_level_indexed_tn5', action='store_true', help='Flag to run assuming that the library is a two-level indexed-TN5 sample.')
	parser.add_argument('--wells_384', action='store_true', help='Flag to run assuming That the known barcode set is the 384 well set.')
	parser.add_argument('--well_ids', action='store_true', help='Flag to output cell IDs that are composed of well IDs rather than the actual sequences.')
	args = parser.parse_args()

	if args.two_level_indexed_tn5 and args.wells_384:
		raise ValueError('There is no 384 well barcode set for indexed Tn5, may not specify both --two_level_indexed_tn5 and --wells_384.')

	indices = args.indices
	name = args.name
	indexsplit = indices.strip().split(':')

	tagi7_indices = indexsplitter(indexsplit[0])
	pcri7_indices = indexsplitter(indexsplit[1])
	pcri5_indices = indexsplitter(indexsplit[2])
	tagi5_indices = indexsplitter(indexsplit[3])

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

	for tagi7_id in tagi7_indices:
		for pcri7_id in pcri7_indices:
			for pcri5_id in pcri5_indices:
				for tagi5_id in tagi5_indices:
					try:
						tagmentation_i7_seq = tagi7[tagi7_id]
						pcr_i7_seq = pcri7[pcri7_id]
						pcr_i5_seq = pcri5[pcri5_id]
						tagmentation_i5_seq = tagi5[tagi5_id]
						
						if args.two_level_indexed_tn5:
							barcodes_string = barcode_to_well.get_two_level_barcode_string(tagmentation_i7_seq, pcr_i7_seq, pcr_i5_seq, tagmentation_i5_seq, bc.nex_two_level_indexed_tn5_to_well, bc.pcr_two_level_indexed_tn5_to_well, args.well_ids)
						else:
							barcodes_string = barcode_to_well.get_barcode_string(tagmentation_i7_seq, pcr_i7_seq, pcr_i5_seq, tagmentation_i5_seq, lig_i7_to_well, lig_i5_to_well, pcr_to_well, args.well_ids)

						print(barcodes_string + '\t' + name)
					except IndexError:
						raise ValueError('One of your specified indices did not fall within the bounds of expected barcodes. If doing a two level indexed tn5 run, this could be because the tn5 barcodes are specified in a 1-8 or 1-12 manner for each end, not 1-96. Or could be some other mistake.')
