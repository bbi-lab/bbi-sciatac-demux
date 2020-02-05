import argparse
import collections
import io_functions as io
import easygrid
import os


# TODO this is duplicated from reanalyze... should consolidate
def load_group_assignments(filename):
    group_assignments = collections.OrderedDict()

    for entry in easygrid.read_delim(filename):
        group_assignments[entry['cell']] = entry['group']

    return group_assignments

if __name__ == '__main__':
    parser = argparse.ArgumentParser('Script to split cells from bam_bed file into multiple bam_bed files according to a grouping sheet.')
    parser.add_argument('bam_bed', help='Input BAM BED file.')
    parser.add_argument('--prefix', required=True, help='Prefix for output files.')
    parser.add_argument('--output_dir', required=True, help='Output directory for files.')
    parser.add_argument('--group_assignments', required=True, help='Tab-delimited file with cell as first column and group as the second column.')
    args = parser.parse_args()

    args = parser.parse_args()
    
    group_assignments = load_group_assignments(args.group_assignments)
    file_handles = {}

    for line in io.open_file(args.bam_bed):
        # Figure out which file should go into
        cell = line.strip().split('\t')[-1]

        assigned_group = group_assignments.get(cell, None)

        if assigned_group is None:
            continue

        # Write to file or open new file as needed
        output_file = os.path.join(args.output_dir, '%s.%s.bed' % (args.prefix, assigned_group))

        if output_file not in file_handles:
            file_handles[output_file] = open(output_file, 'w')

        file_handles[output_file].write(line)

    # Close files
    for filename in file_handles:
        file_handles[filename].close()