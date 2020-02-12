#!/usr/bin/env python3

import sys
import os

nextflowExe = '/net/gs/vol1/home/bge/bin/nextflow'
demuxNextflowScript = '/net/gs/vol1/home/bge/eclipse-workspace/bbi-sciatac-demux/main.nf'
analyzeNextflowScript = '/net/gs/vol1/home/bge/eclipse-workspace/bbi-sciatac-analyze/main.nf'

#
# TODO
#   o  copy samplesheet source file to demux directory
#   o  make a samplesheet file suitable for Andrew's barcode correction utility
#      and write it to the demux directory. Set the samplesheet param in the
#      params file to this converted samplesheet file.
#


valuesInit = { 'procDir'             : { 'index' : 0,
                                         'type' : 'string',
                                         'prompt' : 'Processing (parent) directory',
                                         'optional' : 0,
                                         'default' : '',
                                         'value' : None },
               'runDir'              : { 'index' : 1,
                                         'type' : 'string',
                                         'prompt' : 'Illumina run directory',
                                         'optional' : 0,
                                         'default' : '',
                                         'value' : None },
               'samplesheet'         : { 'index' : 2,
                                         'type' : 'string',
                                         'prompt' : 'Samplesheet file',
                                         'optional' : 0,
                                         'default' : '',
                                         'value' : None },
               'levels'              : { 'index' : 3,
                                         'type' : 'int',
                                         'prompt' : 'Levels',
                                         'optional' : 0,
                                         'default' : 3,
                                         'value' : None },
               'numWells'            : { 'index' : 4,
                                         'type' : 'int',
                                         'prompt' : 'Number of wells',
                                         'optional' : 0,
                                         'default' : 384,
                                         'value' : None },
               'demuxScriptBash'     :   { 'index' : 5,
                                         'type' : 'string',
                                         'prompt' : 'Demux script name',
                                         'optional' : 0,
                                         'default' : 'run.demux.sh',
                                         'value' : None },
               'demuxDir'            : { 'index' : 6,
                                         'type' : 'string',
                                         'prompt' : 'Demux output directory',
                                         'optional' : 0,
                                         'default' : 'demux_dir',
                                         'value' : None },
               'analyzeScriptBash'   : { 'index' : 7,
                                         'type' : 'string',
                                         'prompt' : 'Analyze script name',
                                         'optional' : 0,
                                         'default' : 'run.analyze.sh',
                                         'value' : None },
               'analyzeDir'          : { 'index' : 8,
                                         'type' : 'string',
                                         'prompt' : 'Analyze output directory',
                                         'optional' : 0,
                                         'default' : 'analyze_dir',
                                         'value' : None },
               'paramsFilename'      : { 'index' : 9,
                                         'type' : 'string',
                                         'prompt' : 'Nextflow parameters filename',
                                         'optional' : 0,
                                         'default' : 'params.config',
                                         'value' : None },
               'numCpu'              : { 'index' : 10,
                                         'type' : 'int',
                                         'prompt' : 'Maximum number of CPUs',
                                         'optional' : 0,
                                         'default' : 16,
                                         'value' : None },
               'maxMemoryBcl2fastq'  : { 'index' : 11,
                                         'type' : 'int',
                                         'prompt' : 'Maximum memory for bcl2fastq (GB)',
                                         'optional' : 0,
                                         'default' : 40,
                                         'value' : None },
               'queue'               : { 'index' : 12,
                                         'type' : 'string',
                                         'prompt' : 'Grid engine queue name',
                                         'optional' : 1,
                                         'default' : '',
                                         'value' : None },
               'clusterOptions'      : { 'index' : 13,
                                         'type' : 'string',
                                         'prompt' : 'Grid engine cluster options',
                                         'optional' : 1,
                                         'default' : '',
                                         'value' : None },
               'genomesJson'         : { 'index' : 14,
                                         'type' : 'string',
                                         'prompt' : 'Genomes JSON file',
                                         'optional' : 0,
                                         'default' : '/net/gs/vol1/home/bge/eclipse-workspace/bbi-sciatac-analyze/genomes.json',
                                         'value' : None }
}


def initializeValues( values ):
    for key in values.keys():
        values[key]['value'] = values[key]['default']
    values['procDir']['value'] = os.getcwd()
    return( values )


def getValue( value ):
    if( value['type'] == 'int' ):
        invalue = input( "%s [%d]: " % ( value['prompt'], value['value'] ) )
        invalue = invalue.strip()
        if( invalue != '' ):
            value['value'] = int( invalue )
    elif( value['type'] == 'float' ):
        invalue = input( "%s [%f]: " % ( value['prompt'], value['value'] ) )
        invalue = invalue.strip()
        if( invalue != '' ):
            value['value'] = float( invalue )
    elif( value['type'] == 'string' ):
        invalue = input( "%s [\'%s\']: " % ( value['prompt'], value['value'] ) )
        invalue = invalue.strip()
        if( invalue != '' ):
            value['value'] = invalue
    else:
        print( "unrecognized type value: \'%s\'" % ( value['type'] ), file=sys.stderr )
        sys.exit( -1 )
    return( value )


def getRequiredValues( values ):
    print( '== Enter required values ==' )
    values['procDir'] = getValue( values['procDir'] )
    values['runDir'] = getValue( values['runDir'] )
    values['samplesheet'] = getValue( values['samplesheet'] )
    values['demuxDir']['value'] = values['procDir']['value'] + '/' + values['demuxDir']['value']
    values['analyzeDir']['value'] = values['procDir']['value'] + '/' + values['analyzeDir']['value']
    values['demuxScriptBash']['value'] = values['procDir']['value'] + '/' + values['demuxScriptBash']['value']
    values['analyzeScriptBash']['value'] = values['procDir']['value'] + '/' + values['analyzeScriptBash']['value']
    return( values )


def showValues( values, lvalues, indexFlag ):
    for key in lvalues:
        sindex = '%d) ' % ( values[key]['index'] + 1 ) if indexFlag else ''
        soptional = '(optional) ' if values[key]['optional'] else ''
        if( values[key]['value'] is None ):
            print( '%s%s %s: None' % ( sindex, values[key]['prompt'], soptional ))
        elif( values[key]['type'] == 'int' ):
            print( '%s%s %s: %d' % ( sindex, values[key]['prompt'], soptional, values[key]['value'] ) )
        elif( values[key]['type'] == 'float' ):
            print( '%s%s %s: %f' % ( sindex, values[key]['prompt'], soptional, values[key]['value'] ) )
        elif( values[key]['type'] == 'string' ):
            print( '%s%s %s: \'%s\'' % ( sindex, values[key]['prompt'], soptional, values[key]['value'] ) )
    return( 0 )


def editValues( values, lvalues ):
    print( '== Edit values for run ==' )
    while( 1 ):
        showValues( values, lvalues, 1 )
        print()
        print()
        invalue = input( 'Enter a value index to edit or \'q\' to exit and set up the pipeline: ' )
        invalue = invalue.strip()
        if( invalue == 'q' ):
            return( values )
        else:
            ix = int( invalue ) - 1
            values[lvalues[ix]] = getValue( values[lvalues[ix]] )
            print()
    return( values )


def checkIlluminaRunDirectory( values ):
    isdir = os.path.isdir( values['runDir']['value'] )
    if( not isdir ):
        print( "%s does not exist or is not readable" % ( values['runDir']['value'] ) )
        sys.exit( -1 )
    return( 0 )

    
def makeDirectory( pathName ):
    if( not os.path.exists( pathName ) ):
        os.makedirs( pathName )
    return( 0 )


def writeNextflowConfigFile( values ):
    s = ''
    s += 'process.beforeScript = \". /etc/profile.d/modules.sh\"\n'
    s +='process.executor = \"sge\"\n'
    if( values['queue']['value'] != '' ):
        s += 'process.queue = \"%s\"\n' % ( values['queue']['value'] )
    if( values['clusterOptions']['value'] != '' ):
        s += 'process.clusterOptions = \"%s\"\n' % ( values['clusterOptions']['value'] )
    fn = values['demuxDir']['value'] + '/nextflow.config'
    fp = open( fn, 'w+' )
    print( '%s' % ( s ), file=fp )          
    fp.close()
    fn = values['analyzeDir']['value'] + '/nextflow.config'
    fp = open( fn, 'w+' )
    print( '%s' % ( s ), file=fp )          
    fp.close()
    return( 0 )


def writeNextflowParamsFile( values ):
    s = ''
    s += 'params.run_dir = \"%s\"\n' % ( values['runDir']['value'] )
    s += 'params.demux_dir = \"%s\"\n' % ( values['demuxDir']['value'] )
    s += 'params.analyze_dir = \"%s\"\n' % ( values['analyzeDir']['value'] )
    s += 'params.sample_sheet = \"%s\"\n' % ( values['samplesheet']['value'] )
    s += 'params.genomes_json = \"%s\"\n' % ( values['genomesJson']['value'] )
    s += 'params.num_well = %d\n' % ( values['numWells']['value'] )
    s += 'params.level = %d\n' % ( values['levels']['value'] )
    s += 'params.max_cores = %d\n' % ( values['numCpu']['value'] )
    s += 'params.max_mem_bcl2fastq = %d\n' % ( values['maxMemoryBcl2fastq']['value'] )
    fn = values['demuxDir']['value'] + '/' + values['paramsFilename']['value']
    fp = open( fn, 'w+' )
    print( '%s' % ( s ), file=fp )
    fp.close()
    fn = values['analyzeDir']['value'] + '/' + values['paramsFilename']['value']
    fp = open( fn, 'w+' )
    print( '%s' % ( s ), file=fp )
    fp.close()
    return( 0 )


def writeDemuxScript( values ):
    paramsFilename = values['demuxDir']['value'] + '/' + values['paramsFilename']['value']
    fn = values['demuxScriptBash']['value']
    fp = open( fn, 'w' )
    print( '#!/bin/bash', file=fp )
    print( '', file=fp )
    print( 'NEXTFLOW_EXE=\"%s\"' % ( nextflowExe ), file=fp )
    print( 'NEXTFLOW_SCRIPT=\"%s\"' % ( demuxNextflowScript ), file=fp )
    print( 'PARAMS_FILE=\"%s\"' % ( paramsFilename ), file=fp )
    print( 'WORK_DIR=\"%s\"' % ( values['demuxDir']['value'] + '/work' ), file=fp )
    print( '', file=fp )
    print( 'REPORT_FILE=./dmux.report.html', file=fp )
    print( 'TRACE_FILE=./dmux.trace.tsv', file=fp )
    print( 'TIMELINE_FILE=./dmux.timeline.html', file=fp )
    print( '', file=fp )
    print( 'pushd %s' % ( values['demuxDir']['value'] ), file=fp )
    print( 'date > ./run_start.txt', file=fp )
    print( '', file=fp )
    print( 'PARS="-c $PARAMS_FILE -w $WORK_DIR -with-report $REPORT_FILE -with-trace $TRACE_FILE -with-timeline $TIMELINE_FILE"', file=fp )
    print( '$NEXTFLOW_EXE run $NEXTFLOW_SCRIPT $PARS', file=fp )
    print( '', file=fp )
    print( 'date > ./run_finish.txt', file=fp )
    print( 'popd', file=fp )
    print( '', file=fp )
    fp.close()
    return( 0 )


def writeAnalyzeScript( values ):
    paramsFilename = values['analyzeDir']['value'] + '/' + values['paramsFilename']['value']
    fn = values['analyzeScriptBash']['value']
    fp = open( fn, 'w' )
    print( '#!/bin/bash', file=fp )
    print( '', file=fp )
    print( 'NEXTFLOW_EXE=\"%s\"' % ( nextflowExe ), file=fp )
    print( 'NEXTFLOW_SCRIPT=\"%s\"' % ( analyzeNextflowScript ), file=fp )
    print( 'PARAMS_FILE=\"%s\"' % ( paramsFilename ), file=fp )
    print( 'WORK_DIR=\"%s\"' % ( values['analyzeDir']['value'] + '/work' ), file=fp )
    print( '', file=fp )
    print( 'REPORT_FILE=./analyze.report.html', file=fp )
    print( 'TRACE_FILE=./analyze.trace.tsv', file=fp )
    print( 'TIMELINE_FILE=./analyze.timeline.html', file=fp )
    print( '', file=fp )
    print( 'pushd %s' % ( values['analyzeDir']['value'] ), file=fp )
    print( 'date > ./run_start.txt', file=fp )
    print( '', file=fp )
    print( 'PARS="-c $PARAMS_FILE -w $WORK_DIR -with-report $REPORT_FILE -with-trace $TRACE_FILE -with-timeline $TIMELINE_FILE"', file=fp )
    print( '$NEXTFLOW_EXE run $NEXTFLOW_SCRIPT $PARS', file=fp )
    print( '', file=fp )
    print( 'date > ./run_finish.txt', file=fp )
    print( 'popd', file=fp )
    print( '', file=fp )
    fp.close()
    return( 0 )


# initialize values from the defaults
values = valuesInit
values = initializeValues( values )

# gather some supporting information
nvalues = len( values.keys() )
lvalues = list( range( nvalues ) )
for key in values.keys():
    lvalues[values[key]['index']] = key

# get required run-specific information
values = getRequiredValues( values )
print()
print()

# allow the user to edit values
values = editValues( values, lvalues )

# check for Illumina run directory.
checkIlluminaRunDirectory( values )

# make demux and analysis directories.
makeDirectory( values['demuxDir']['value'] )
makeDirectory( values['analyzeDir']['value'] )

# write Nextflow configuration files.
writeNextflowConfigFile( values )

# write Nextflow parameter files.
writeNextflowParamsFile( values )

# write samplesheet file
# I need more information in order to make the samplesheet file...
#writeSamplesheetFile( values )

# write demux script.
writeDemuxScript( values )

# write analysis script.
writeAnalyzeScript( values )

