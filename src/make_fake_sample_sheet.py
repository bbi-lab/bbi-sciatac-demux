#!/usr/bin/env python
# Andrew's sample sheet creator

import os
import argparse

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))


def get_sample_sheet_text(p7_index_length, p5_index_length):
    """
    Gets the sample sheet text that will demux cells into one set of files
    """

    sample_sheet_template = """[DATA]
Lane,Sample_ID,Sample_Name,index,index2
%s"""

    line = ',fake,fake,%s,%s' % ('N' * p7_index_length, 'N' * p5_index_length)
    return sample_sheet_template % line


if __name__ == '__main__':
    # Script
    parser = argparse.ArgumentParser('Script make sample sheet')

    # Required args
    parser.add_argument('--p7_index_length', required=True, help='Length of p7 index.')
    parser.add_argument('--p5_index_length', required=True, help='Length of p5 index.')
    args = parser.parse_args()

    # Set up samplesheet for BCL2FASTQ
    sample_sheet_text = get_sample_sheet_text( int(args.p7_index_length), int(args.p5_index_length) )
    sample_sheet_path = os.path.join('SampleSheet.csv')
    with open(sample_sheet_path, 'w') as sample_sheet:
        sample_sheet.write(sample_sheet_text)


