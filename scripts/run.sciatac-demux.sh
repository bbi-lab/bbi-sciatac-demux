#!/bin/bash

#
# Edit this file to set the NEXTFLOW variable to the location
# of the Nextflow executable, and NF_MAIN to the location of
# the bbi-sciatac-demux Nextflow pipeline script, main.nf.
#

#
# Copy this file to the directory where you will run the pipeline.
# The experiment.config file must be in the directory where you
# run this script.
#

#
# Nextflow executable and pipeline script locations.
# Note: set the paths in the three variables below.
#
NEXTFLOW="$HOME/bin/nextflow"
NF_HOME="$HOME/git/bbi-sciatac-demux"
NF_MAIN="${NF_HOME}/main.nf"

#
# Current date and time.
#
NOW=`date '+%Y%m%d_%H%M%S'`

#
# Path to the Nextflow processing run configuration file.
#
PWD=`pwd`
CONFIG_FILE="$PWD/experiment.config"

#
# Get the path to the demux output directory from
# the configuration file and set the Nextflow work
# directory to be in the demux output directory.
# Notes:
#   o  bbi-sciatac-demux/main.nf and bbi-scripts-analyze/main.nf also define
#      DEMUX_DIR as ${params.output_dir}/demux_out so changing DEMUX_DIR here
#      requires changing those files as well.
#
OUTPUT_DIR=`cat $CONFIG_FILE | sed 's/[ ]*//g' | awk 'BEGIN{FS="="}{if($1=="params.output_dir"){print$2}}' | sed 's/"//g'`
DEMUX_DIR="$OUTPUT_DIR/demux_out"
WORK_DIR="$DEMUX_DIR/work"

#
# Use REPORT_FIL, TRACE_FIL, and TIMELINE_FIL so that
# Nextflow writes some informative processing information
# in the analyze output directory.
#
# REPORT_FIL=$DEMUX_DIR/run_reports/demux.report.${NOW}.html
TRACE_FIL=$DEMUX_DIR/run_reports/demux.trace.${NOW}.tsv
# TIMELINE_FIL=$DEMUX_DIR/run_reports/demux.timeline.${NOW}.html

#
# Nextflow run parameters.
#
# PARS="-c $CONFIG_FILE -w $WORK_DIR -with-report $REPORT_FIL -with-trace $TRACE_FIL -with-timeline $TIMELINE_FIL -resume"
PARS="-c $CONFIG_FILE -w $WORK_DIR -with-trace $TRACE_FIL -resume"

mkdir -p $DEMUX_DIR/run_reports
pushd $DEMUX_DIR

date > run_reports/run_start.${NOW}.txt

#
# Run Nextflow sci-ATAC demux pipeline.
#
$NEXTFLOW run $NF_MAIN $PARS

date > run_reports/run_finish.${NOW}.txt

popd

#
# Set run directory file and directory permissions.
#
${NF_HOME}/scripts/set_run_permissions.sh ${OUTPUT_DIR}

