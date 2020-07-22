#!/usr/bin/python3

import os
import xml.etree.ElementTree as ET
import glob
import sys
import gzip
import re

NEXTSEQ = 'NextSeq'
MISEQ = 'MiSeq'
NOVASEQ = 'NovaSeq'
HISEQ4000 = 'HiSeq4000'
HISEQ3000 = 'HiSeq3000'
HISEQ = 'HiSeq'
UNKNOWN_SEQUENCER = 'unknown'

SEQUENCERS_P5_RC_MAP = {
    NEXTSEQ: True,
    MISEQ: False,
    NOVASEQ: False,
    HISEQ4000: True,
    HISEQ3000: False
}

def open_file(f, mode='rt'):
    if f.endswith('.gz'):
        return gzip.open(f, mode)
    else:
        return open(f, mode)


# Notes on MiSeq runParameters.xml
#   o  a MiSeq runParameters.xml file differs in comparison to a NextSeq RunParameters.xml
#      file in at least the following ways
#        o  file name differs as written above
#        o  FlowcellRFIDTag vs FlowCellRfidTag
#        o  ScannerID vs InstrumentID
#        o  Read information format differs
#             o  MiSeq: in <Reads> block
#             o  NextSeq: in <Setup> block
#
def get_run_info(flow_cell_path):
    """
    Helper function to get some info about the sequencing runs.
    Args:
        flow_cell_path (str): Path to BCL directory for run.
    Returns:
        dict: basic statistics about run, like date, instrument, number of lanes, flowcell ID, read lengths, etc.
    """
    
    print( 'get_run_info: pathname: %s' % ( flow_cell_path ),file=sys.stderr)
    
    run_stats = {}

    bcl_run_info = os.path.join(flow_cell_path, 'RunParameters.xml*')
    print( 'bcl_run_info1: %s' % ( bcl_run_info ),file=sys.stderr)
    bcl_run_info = glob.glob(bcl_run_info)
    print( 'bcl_run_info2: %s' % ( bcl_run_info ),file=sys.stderr)
    if not bcl_run_info:
        raise ValueError('BCL RunParameters.xml not found for specified flowcell: %s' % bcl_run_info)
    else:
        bcl_run_info = bcl_run_info[0]

    # Set up a few nodes for parsing
    tree = ET.parse(open_file(bcl_run_info))

    setup_node = tree.getroot().find("Setup")
    if setup_node is None:
        setup_node = tree.getroot()

    flowcell_node = tree.getroot().find("FlowCellRfidTag")
    instrument_id_node = tree.getroot().find('InstrumentID')
    run_start_date_node = tree.getroot().find('RunStartDate')

    # Now actually populate various stats
    run_stats['flow_cell_id'] = flowcell_node.find('SerialNumber').text
    run_stats['date'] = run_start_date_node.text
    run_stats['instrument'] = instrument_id_node.text
    run_stats['lanes'] = int(setup_node.find('NumLanes').text)
    
    run_stats['r1_length'] = int(setup_node.find('Read1').text)
    run_stats['r2_length'] = int(setup_node.find('Read2').text)
    run_stats['p7_index_length'] = int(setup_node.find('Index1Read').text)
    run_stats['p5_index_length'] = int(setup_node.find('Index2Read').text)

    application = setup_node.find('ApplicationName').text
    application_version = setup_node.find('ApplicationVersion')
    if NEXTSEQ in application:
        run_stats['instrument_type'] = NEXTSEQ
    elif MISEQ in application:
        run_stats['instrument_type'] = MISEQ
    elif NOVASEQ in application:
        run_stats['instrument_type'] = NOVASEQ
    elif HISEQ in application:
        app_string = re.search(r'[\d\.]+', application_version).group()
        app_major_version = int(app_string.split('.')[0])

        if app_major_version > 2:
            run_stats['instrument_type'] = HISEQ4000
        else:
            run_stats['instrument_type'] = HISEQ3000
    else:
        run_stats['instrument_type'] = UNKNOWN_SEQUENCER

    return run_stats


def reverse_complement_i5(name):
    """
    Take a BCL directory or instrument type (NextSeq, MiSeq, NovaSeq, HiSeq4000, HiSeq3000) and return whether or not i5 should be reverse complemented.
    This assumes that NextSeq instruments and other similar machines should be reverse complemeted whereas MiSeq should not.
    Args:
        name (str): BCL directory or one of the instrument types as mentioned above    
    
    Returns:
        bool: True if user would typically reverse complement i5 index and False otherwise.
    """
    
    if name in SEQUENCERS_P5_RC_MAP:
        sequencer_type = name
    elif os.path.exists(name):
        sequencer_type = get_run_info(name)['instrument_type']
        
        if sequencer_type not in SEQUENCERS_P5_RC_MAP:
            raise ValueError('Sequencer type detected from BCL is %s, which is not in our known list of which sequencers require P5 reverse complementing or not.' % sequencer_type)
    else:
        raise ValueError('Invalid input, could not detect BCL or instrument ID.')

    return SEQUENCERS_P5_RC_MAP[sequencer_type]

if __name__ == '__main__':
    pathname = sys.argv[1]
#    pathname = '/net/shendure/vol9/seq/NEXTSEQ/190801_NS500488_0863_AHN3M7AFXY'
    run_stats = get_run_info(pathname)
    
    for key in run_stats.keys():
        print( '%s %s' % ( key, run_stats[key] ), file=sys.stdout )
    print( 'reverse_complement_i5 %d' % ( int( reverse_complement_i5( run_stats['instrument_type'] ) == True ) ) )
    