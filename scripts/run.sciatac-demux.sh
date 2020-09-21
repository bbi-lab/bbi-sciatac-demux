#!/bin/bash

#
# Path to the Nextflow processing run configuration file.
#
CONFIG_FILE="/xxx/params.config"

#
# Nextflow executable and pipeline script locations.
#
NEXTFLOW="/xxx/bin/nextflow"
NF_DEMUX="/xxx/bbi-sciatac-demux/main.nf"

#
# Get the path to the demux output directory from
# the configuration file and set the Nextflow work
# directory to be in the demux output directory.
#
DEMUX_DIR=`cat $CONFIG_FILE | sed 's/[ ]*//g' | awk 'BEGIN{FS="="}{if($1=="params.demux_dir"){print$2}}' | sed 's/"//g'`
WORK_DIR="$DEMUX_DIR/work"

#
# Use REPORT_FIL, TRACE_FIL, and TIMELINE_FIL so that
# Nextflow writes some informative processing information
# in the demux output directory.
#
REPORT_FIL=$DEMUX_DIR/demux.report.html
TRACE_FIL=$DEMUX_DIR/demux.trace.tsv
TIMELINE_FIL=$DEMUX_DIR/demux.timeline.html

#
# Nextflow run parameters.
# Note: I include -resume for convenience. I does not affect the initial run.
#
PARS="-c $CONFIG_FILE -w $WORK_DIR -with-report $REPORT_FIL -with-trace $TRACE_FIL -with-timeline $TIMELINE_FIL -resume"

mkdir -p $DEMUX_DIR
pushd $DEMUX_DIR

date > ./run_start.txt

#
# Run Nextflow sci-ATAC demux pipeline.
#
$NEXTFLOW run $NF_DEMUX $PARS

date > run_finish.txt

popd
