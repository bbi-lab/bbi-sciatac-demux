/*
** This pipeline is written for Nextflow DSL 1.
*/
nextflow.enable.dsl = 1

/* Default configuration parameters.
** Typically, many of these are given in a configuration file.
** The configuration file is specified as a parameter on the
** NextFlow command, for example,
**   nextflow run main.nf -c nextflow_parameters.config
**
** Run in grid qlogin shell with at least 10G of memory (Andrew's allocation for barcode correction)
** Notes:
**   o  save intermediate files option?
**   o  think about how to deal with python and pypy virtual environment
**   o  consider writing a script that builds files required for the run
**      and sets up the output/run directory. Files created for run
**
**        o  samplesheet
**        o  nextflow parameters file
**        o  nextflow configuration file
**        o  bash script to run this Nextflow script
**      
** Suggested command line parameters to use when running this script
**   nextflow run main.nf -w <work_dirname> -with-report <report_filename> -with-trace <trace_filename> -with-timeline <timeline_filename>
**
** Notes:
**   o  see bbi-sciatac-analyze/main.nf header comments for notes
**      aspects of Groovy and Nextflow. For example, variable
**      scope in Groovy.
*/


import java.nio.file.Files
import java.nio.file.Path 
import java.nio.file.Paths
import groovy.json.JsonOutput
import groovy.json.JsonSlurper


/*
** Nextflow and main.nf versions.
** manifest.version in set in the nextflow.config file.
*/
nextflow_version = nextflow.version.toString()
bbi_sciatac_demux_version = "$params.version"


/*
** Run date/time.
*/
def timeNow = new Date()


/*
** Minimum required samplesheet file version.
*/
minimum_samplesheet_version = "3.0.0"


/*
** Where to find scripts.
** Note: script_dir needs to be visible within Groovy functions
**       so there is no 'def', which makes it global.
*/
pipeline_path="$workflow.projectDir"
script_dir="${pipeline_path}/src"


/*
** Set errorStrategy directive policy.
**
** Allowed values:
**   "terminate"  (default NextFlow policy)
**   "finish"
**   "ignore"
**   "retry"      (if set to retry, NextFlow retries 1 time (by default: can be changed with 'maxRetries' directive)
**
** Notes:
**   o  the NextFlow documentation describes dynamic retries: errorStrategy { sleep(Math.pow(2, task.attempt) * 200); return 'retry' }
**      this increases exponentially the time between retries, in case of long latencies
*/
def onError = { return( "retry" ) }

/*
** ================================================================================
** Initial set up and checks.
** ================================================================================
*/

/*
** Initial pre-defined, required parameter values.
** Notes:
**   o  params.index_recipe uses the implicit index recipes when set to 0.
**      Use the default unless you know what you are doing.
*/
params.help = false
params.bcl2fastq_cpus = 6
params.max_mem_bcl2fastq = 40
params.demux_buffer_blocks = 8192
params.index_recipe = 0

/*
** Initialize optional parameters to null.
*/
// none at this time

/*
** Define/initialize some internal parameters.
** Notes:
*/
def sample_sheet = params.sample_sheet
def run_dir      = params.run_dir
def demux_buffer_blocks = params.demux_buffer_blocks

def options_barcode_correct = ''

/*
** Print usage when run with --help parameter.
*/
if (params.help) {
    writeHelp()
    exit( 0 )
}

/*
** Check for required parameters.
*/
if( !params.run_dir ) {
	printErr( "Error: missing params.run_dir in parameters file" )
	System.exit( -1 )
}
if( !params.output_dir ) {
	printErr( "Error: missing params.output_dir in parameters file" )
	System.exit( -1 )
}
if( !params.sample_sheet ) {
	printErr( "Error: missing params.sample_sheet in parameters file" )
	System.exit( -1 )
}

/*
** Global variables accessible in process blocks.
*/
output_dir = params.output_dir.replaceAll("/\\z", "")
demux_dir = output_dir + '/demux_out'
log_dir = output_dir + '/demux_log_dir'
tmp_dir = output_dir + '/tmp'

/*
** Check sample sheet file.
*/
checkSamplesheet( params )

/*
** Check that required directories exist or can be made.
*/
checkDirectories( params, log_dir, tmp_dir )

/*
** Read sample sheet file.
*/
sampleSheetMap = readSampleSheetJson( params )

/*
** Test samplesheet version.
*/
if( ! checkSamplesheetVersion( sampleSheetMap['json_file_version'], minimum_samplesheet_version ) ) {
  println "Error: bad samplesheet version"
  System.exit( -1 )
}

/*
** Report run parameter values.
*/
reportRunParams( params, sampleSheetMap )

/*
** Archive configuration and samplesheet files in demux_dir.
*/
archiveRunFiles( params, timeNow )

/*
** Make a groovy list of sample names.
*/
sampleList = makeSampleList(sampleSheetMap)

/*
** Save workflow.runName to a file that can be
** read by the logger in each process.
** Note: using ${workflow.runName} in a process
** causes the -resume to fail because the runName
** changes from run-to-run.
*/
File tfile = new File("${tmp_dir}/nextflow_run_name.txt")
tfile.write("${workflow.runName}")

/*
** Read Illumina run information.
*/
illuminaRunInfoMap = readIlluminaRunInfo( params )

/*
** Check that these are paired-end reads.
*/
if( illuminaRunInfoMap['paired_end'] == false )
{
  throw new Exception('Single-end reads detected: paired-end reads required')
}

/*
** Write run information to args.json file.
*/
writeArgsJson( params, timeNow )

/*
** Does the i5 sequence require reverse complementing?
** If yes, set sequence_flag.
*/
def sequencer_flag = ( illuminaRunInfoMap['reverse_complement_i5'] ) ? '-X' : ''

/*
** Add optional barcode correction parameters.
*/
if( sampleSheetMap['tn5_barcodes'] && ( sampleSheetMap['level'] != 2 || sampleSheetMap['number_wells'] != 96 ) ) {
    throw new Exception('tn5_barcode requires level == 2 and number_wells == 96')
}
if( sampleSheetMap['level'] == 2 && sampleSheetMap['number_wells'] == 96 && sampleSheetMap['tn5_barcodes'] ) {
	options_barcode_correct += ' --two_level_indexed_tn5'
}
if( sampleSheetMap['number_wells'] == 384 ) {
	options_barcode_correct += ' --wells_384'
}
if( sampleSheetMap['use_all_barcodes'] ) {
  options_barcode_correct += ' --no_mask'
}
if( params.index_recipe > 0 ) {
  options_barcode_correct += sprintf("  --index_recipe %d", params.index_recipe)
}

/*
** Use well ids as read names. This is required for downstream
** quality control evaluation.
*/
options_barcode_correct += ' --well_ids'

/*
** Set number of threads and amount of memory per core
** allocated for bcl2fastq run.
*/
if( params.bcl2fastq_cpus > 16 ) {
	num_threads_bcl2fasta_process = 16
	mem_bcl2fastq = params.max_mem_bcl2fastq / 16
} else {
	num_threads_bcl2fasta_process = params.bcl2fastq_cpus
	mem_bcl2fastq = params.max_mem_bcl2fastq / num_threads_bcl2fasta_process
}

if( num_threads_bcl2fasta_process / 2 < 4 ) {
	num_threads_bcl2fastq_io = num_threads_bcl2fasta_process / 2
} else {
	num_threads_bcl2fastq_io = 4
}


/*
** Write nextflow and main.nf versions to a log file.
*/
process log_pipeline_versions {

    script:
    """
    # bash watch for errors
    set -ueo pipefail

    PROCESS_BLOCK='log_pipeline_versions'
    SAMPLE_NAME="pipeline"
    START_TIME=`date '+%Y%m%d:%H%M%S'`
    STOP_TIME=`date '+%Y%m%d:%H%M%S'`
    PIPELINE_VERSIONS_JSON="{\\\"nextflow_version\\\": ${nextflow_version}, {\\\"bbi_sciatac_demux_version\\\": ${bbi_sciatac_demux_version}}"

    $script_dir/pipeline_logger.py \
    -r `cat ${tmp_dir}/nextflow_run_name.txt` \
    -n \${SAMPLE_NAME} \
    -p \${PROCESS_BLOCK} \
    -v "echo -n \\\"Nextflow version: ${nextflow_version}\\\"" "echo -n \\\"bbi_sciatac_demux_version: ${bbi_sciatac_demux_version}\\\"" \
    -s \${START_TIME} \
    -e \${STOP_TIME} \
    -d ${log_dir}
    """
}


/*
** Make a fake sample sheet required by bcl2fastq. The sample sheet has Ns for
** the sequence in order to force bcl2fastq to make undetermined fastqs.
*/
process make_fake_sample_sheet {
  cache 'lenient'
  errorStrategy onError

  output:
  file("SampleSheet.csv") into makeFakeSampleSheetOutChannel

  script:
  """
  # bash watch for errors
  set -ueo pipefail

  $script_dir/make_fake_sample_sheet.py --p7_index_length=${illuminaRunInfoMap['p7_index_length']} --p5_index_length=${illuminaRunInfoMap['p5_index_length']}
  """
}


/*
** Run bcl2fastq to make fastq files from Illumina bcl files.
**
** Notes:
**   o  copy 'Reports' and 'Stats' directories to an accessible directory.
*/
makeFakeSampleSheetOutChannel
    .set { bcl2fastqInChannel }

process bcl2fastq {
  cache 'lenient'
  errorStrategy onError
  cpus num_threads_bcl2fasta_process
  memory "$mem_bcl2fastq GB"    
//  note: the following publishDir statement may be commented out in bbi_dmux. I want it to work for diagnostics during development.
  publishDir path: "$demux_dir/fastqs_bcl2fastq", pattern: "Undetermined_S0_*.fastq.gz", mode: 'copy'
  publishDir path: "$demux_dir/fastqs_bcl2fastq", pattern: "Stats/*", mode: 'copy'
  publishDir path: "$demux_dir/fastqs_bcl2fastq", pattern: "Reports/*", mode: 'copy'

  input:
    file(inFile) from bcl2fastqInChannel
//    file(inFile) from makeFakeSampleSheetOutChannel

  output:
    set file("Undetermined_S0_*_R1_001.fastq.gz"), file("Undetermined_S0_*_R2_001.fastq.gz") into bcl2fastq_fastqsOutChannel mode flatten
    file("Stats/*") into  bcl2fastq_stats
    file("Reports/*") into bcl2fastq_reports

    /*
    ** Notes from bcl2fastq manual:
    **   o  sciatac_pipeline sets interop-dir to the same path as output-dir
    **   o  bcl2fastq defaults
    **        o  --loading-threads: 4
    **        o  --writing-threads: 4 (must be <= # samples. If you specify more threads than samples, the
    **                                 extra threads do no work but cost time due to context switching.)
    **        o  --processing-threads: all available cores/threads
    */
    /*
    ** Notes on development (bge)
    **   o  add --tiles option to restrict conversion during development
    **        To select all tiles ending with 5 in all lanes:
    **        --tiles [0–9][0–9][0–9]5
    */
    script:
    """
    # bash watch for errors
    set -ueo pipefail

    PROCESS_BLOCK='bcl2fastq'
    SAMPLE_NAME="lane"
    START_TIME=`date '+%Y%m%d:%H%M%S'`

    # bcl2fastq command line parameters
    #   ignore-missing-positions: The software ignores missing or corrupt cluster location
    #                             files. When cluster location files are missing, the
    #                             software writes unique coordinate positions into the FASTQ
    #                             header.
    #   ignore-missing-filter: The software ignores missing or corrupt filter files and assumes
    #                          that all clusters in tiles with missing filter files passed filter.
    #   ignore-missing-bcls: The software ignores missing or corrupt BCL files and assumes 'N'/'#'
    #                        for missing calls.
    #   ignore-missing-controls: Missing or corrupt control files are ignored. Missing
    #                            controls: 0
    #
    bcl2fastq --runfolder-dir      ${params.run_dir} \
              --output-dir         . \
              --interop-dir        . \
              --sample-sheet       $inFile \
              --loading-threads    ${num_threads_bcl2fastq_io} \
              --processing-threads ${num_threads_bcl2fasta_process}  \
              --writing-threads    ${num_threads_bcl2fastq_io} \
              --ignore-missing-positions \
              --ignore-missing-controls \
              --ignore-missing-filter \
              --ignore-missing-bcls

    $script_dir/json_extractor.py \
    -i Stats/Stats.json \
    -k ReadInfosForLanes ConversionResults \
    -o bcl2fastq_stats.json

    STOP_TIME=`date '+%Y%m%d:%H%M%S'`
    $script_dir/pipeline_logger.py \
    -r `cat ${tmp_dir}/nextflow_run_name.txt` \
    -n \${SAMPLE_NAME} \
    -p \${PROCESS_BLOCK} \
    -v 'bcl2fastq --version' \
    -s \${START_TIME} \
    -e \${STOP_TIME} \
    -J bcl2fastq_stats.json \
    -d ${log_dir} \
    -c "bcl2fastq --runfolder-dir      ${params.run_dir} \
--output-dir         . \
--interop-dir        . \
--sample-sheet       $inFile \
--loading-threads    ${num_threads_bcl2fastq_io} \
--processing-threads ${num_threads_bcl2fasta_process}  \
--writing-threads    ${num_threads_bcl2fastq_io} \
--ignore-missing-positions \
--ignore-missing-controls \
--ignore-missing-filter \
--ignore-missing-bcls"
    """
    //               --tiles s_1 \
}

bcl2fastq_fastqsOutChannel
   .into { bcl2fastq_fastqsOutChannelCopy01;
           bcl2fastq_fastqsOutChannelCopy02 }


/*
** Run bar code correction script (barcode_correct_sciatac.py).
**
** Notes:
**   o  need to activate pypy environment and then deactivate when done
**   o  copy 'out.*.correction_stats.json' to an accessible directory
**   o  I suspect that splitting fastq files and combining them uses more
**      time than one saves by distributing the barcode correction. I base
**      this on the observation that the file compressions takes about the
**      same amount of time as the barcode correction. I assume that the
**      two scale similarly with file size.
**   o  barcode_correct_sciatac.py uses pigz to compress fastq files. pigz is
**      supposed to use all available processors to compress files.
**   o  the barcode_stats_json output channel is a dummy channel. It appears
**      that the files do not get copied by the publishDir
**      directive if these files are not in a channel.
*/
process barcode_correct {
  cache 'lenient'
  errorStrategy onError
  publishDir    path: "${demux_dir}", saveAs: { qualifyFilename( it, "fastqs_barcode" ) }, pattern: "*.fastq.gz", mode: 'copy'
  publishDir    path: "${demux_dir}/fastqs_barcode", pattern: "*.stats.json", mode: 'copy'
  publishDir    path: "${demux_dir}/fastqs_barcode", pattern: "*.index_counts.csv", mode: 'copy'
  publishDir    path: "${demux_dir}/fastqs_barcode", pattern: "*.tag_pair_counts.csv", mode: 'copy'
  publishDir    path: "${demux_dir}/fastqs_barcode", pattern: "*.pcr_pair_counts.csv", mode: 'copy'
  
  input:
    // set file(R1), file(R2) from bcl2fastq_fastqsOutChannelCopy01.splitFastq(by: params.fastq_chunk_size, file: true, pe: true)
    set file(R1), file(R2) from bcl2fastq_fastqsOutChannelCopy01
    
  output:
    file "*.fastq.gz" into barcode_fastqs mode flatten
    file "*.stats.json" into barcode_stats_json mode flatten
    file "*.index_counts.csv" into index_counts_csv mode flatten
    file "*.tag_pair_counts.csv" into tag_pair_counts_csv mode flatten
    file "*.pcr_pair_counts.csv" into pcr_pair_counts_csv mode flatten
    file ("stop_flag") into bcl2fastq_fastqsStopFlagOutChannel

  script:
  """
  # bash watch for errors
  set -ueo pipefail

  PROCESS_BLOCK='barcode_correct'
  SAMPLE_NAME="lane"
  START_TIME=`date '+%Y%m%d:%H%M%S'`
  LANE_ID=`echo ${R1} | awk 'BEGIN{FS="_"}{print\$3}'`

  source $pipeline_path/load_pypy_env_reqs.sh
  PS1=\${PS1:-}
  source $script_dir/pypy_env/bin/activate

  #
  # Partition reads into sample subsets.
  #
  pypy $script_dir/barcode_correct_sciatac.py \
                       --samplesheet $sample_sheet \
                       -1 <(zcat $R1) \
                       -2 <(zcat $R2) \
                       --filename $R1 \
                       --out_dir . \
                       --stats_out 1 \
                       --num_pigz_threads ${task.ext.num_pigz_threads} \
                       --write_buffer_blocks ${demux_buffer_blocks} \
                       $options_barcode_correct \
                       $sequencer_flag
  deactivate

  STOP_TIME=`date '+%Y%m%d:%H%M%S'`
  $script_dir/pipeline_logger.py \
  -r `cat ${tmp_dir}/nextflow_run_name.txt` \
  -n  \${SAMPLE_NAME} \
  -x \${LANE_ID} \
  -p \${PROCESS_BLOCK} \
  -v 'echo "not available"' \
  -s \${START_TIME} \
  -e \${STOP_TIME} \
  -f *.stats.json \
  -d ${log_dir} \
  -c "pypy $script_dir/barcode_correct_sciatac.py \
--samplesheet $sample_sheet \
-1 <(zcat $R1) \
-2 <(zcat $R2) \
--filename $R1 \
--out_dir . \
--stats_out 1 \
--num_pigz_threads ${task.ext.num_pigz_threads} \
--write_buffer_blocks ${demux_buffer_blocks} \
$options_barcode_correct \
$sequencer_flag"

  touch stop_flag
  """
}


/*
** This is the end of the 'fastq-producing' sciatac_pipeline processing steps.
** However, this script requires still statistics/report generating functions.
*/

/*
** Run fastqc on bcl2fastq output (lanes) files.
*/
process fastqc_lanes {
  cache 'lenient'
  errorStrategy onError
  publishDir path: "$output_dir", pattern: "fastqc_lanes", mode: 'copy'

  input:
    file fastq from bcl2fastq_fastqsOutChannelCopy02.collect()

  output:
    file fastqc_lanes into fastqcLanesOutChannel
    file ("stop_flag") into fastqcLanesStopFlagOutChannel

  script:
  """
  # bash watch for errors
  set -ueo pipefail

  PROCESS_BLOCK='fastqc_lanes'
  SAMPLE_NAME="lane"
  START_TIME=`date '+%Y%m%d:%H%M%S'`

  LANE_ID=`echo ${fastq} | awk 'BEGIN{FS="_"}{print\$3}'`

  mkdir fastqc_lanes
  fastqc *.fastq.gz -t $task.cpus -o fastqc_lanes

  STOP_TIME=`date '+%Y%m%d:%H%M%S'`

  $script_dir/pipeline_logger.py \
  -r `cat ${tmp_dir}/nextflow_run_name.txt` \
  -n \${SAMPLE_NAME} \
  -x \${LANE_ID} \
  -p \${PROCESS_BLOCK} \
  -v 'fastqc -version' \
  -s \${START_TIME} \
  -e \${STOP_TIME} \
  -d ${log_dir} \
  -c "fastqc *.fastq.gz -t $task.cpus -o fastqc_lanes"

  touch stop_flag
  """
}


/*
** Run fastqc on barcode corrected fastq files by sample.
*/
process fastqc_samples {
  cache 'lenient'
  errorStrategy onError
  publishDir path: "$output_dir", pattern: "fastqc_sample", mode: 'copy'

  input:
    file fastq from barcode_fastqs.collect()

  output:
    file fastqc_sample into fastqcSampleOutChannel
    file ("stop_flag") into fastqcSampleStopFlagOutChannel

  script:
  """
  # bash watch for errors
  set -ueo pipefail

  sample_name=`echo "${fastq}" | awk 'BEGIN{FS="-"}{print\$1}'`

  PROCESS_BLOCK='fastqc_samples'
  SAMPLE_NAME="all_samples"
  START_TIME=`date '+%Y%m%d:%H%M%S'`


  mkdir fastqc_sample
  fastqc *.fastq.gz -t $task.cpus -o fastqc_sample

  STOP_TIME=`date '+%Y%m%d:%H%M%S'`
  $script_dir/pipeline_logger.py \
  -r `cat ${tmp_dir}/nextflow_run_name.txt` \
  -n \${SAMPLE_NAME} \
  -p \${PROCESS_BLOCK} \
  -v 'fastqc -version' \
  -s \${START_TIME} \
  -e \${STOP_TIME} \
  -d ${log_dir} \
  -c "fastqc *.fastq.gz -t $task.cpus -o fastqc_sample"

  touch stop_flag
  """
}


/*
** Make demux dash files.
**   Required input files:
**     o  args.json file: this is made outside the pipeline, initially
**     o  files '%s/fastqs_barcode/RUN001_L%03d.stats.json' typically in <run_process_dir>/demux_out/fastqs_barcode directory
**     o  files '/RUN001_', lane, '.pcr_pair_counts.csv', typically in <run_process_dir>/demux_out/fastqs_barcode directory
**     o  files '/RUN001_', lane, '.index_counts.csv', typically in <run_process_dir>/demux_out/fastqs_barcode directory
**     o  process barcode_correct output channels
**          o  file "*.stats.json" into barcode_stats_json mode flatten
**          o  file "*.index_counts.csv" into index_counts_csv mode flatten
**          o  file "*.tag_pair_counts.csv" into tag_pair_counts_csv mode flatten
**          o  file "*.pcr_pair_counts.csv" into pcr_pair_counts_csv mode flatten
**   Runs scripts:
**     o  make_dashboard_data.R (Rscript)
**     o  make_run_data.py
**   Requires additional files/directories:
**     o  demux_dash/js/demux.js and footer.js and style and ...
**   Writes to
**     o  demux_dash/js/run_data.js
**     o  demux_dash/js/recover_summary.js  (eventually)
**     o  demux_dash/img/*
**
**   Input channels:
**     process: barcode_correct
**       output:
**         file "*.stats.json" into barcode_stats_json mode flatten
**         file "*.index_counts.csv" into index_counts_csv mode flatten
**         file "*.tag_pair_counts.csv" into tag_pair_counts_csv mode flatten
**         file "*.pcr_pair_counts.csv" into pcr_pair_counts_csv mode flatten
**
**   Output channels:
**     process demux_dash (!)
**
**   Questions:
**     o  how do I transfer all files in input channel to one process?
*/
project_name = output_dir.substring(output_dir.lastIndexOf("/")+1);

process demux_dash {
  cache 'lenient'
  errorStrategy onError
  publishDir path: "$output_dir", pattern: "demux_dash", mode: 'copy'

  input:
    file stats_json from barcode_stats_json.collect()
    file index_counts from index_counts_csv.collect()
    file tag_pair_counts from tag_pair_counts_csv.collect()
    file pcr_pair_counts from pcr_pair_counts_csv.collect()

  output:
    file demux_dash into demux_dashOutChannel
    file ("stop_flag") into demux_dashStopFlagOutChannel

  script:
  """
  # bash watch for errors
  set -ueo pipefail

  PROCESS_BLOCK='demux_dash'
  SAMPLE_NAME="dashboard"
  START_TIME=`date '+%Y%m%d:%H%M%S'`

  mkdir demux_dash
  cp -R $script_dir/skeleton_dash/* demux_dash

  $script_dir/make_dashboard_data.R  --input_folder="." \
                                     --samplesheet=$sample_sheet \
                                     --project_name=$project_name \
                                     --image_output_folder="demux_dash/img"

  $script_dir/make_run_data.py --input_folder="." \
                               --input_file="$demux_dir/args.json" \
                               --output_file="demux_dash/js/run_data.js"

  STOP_TIME=`date '+%Y%m%d:%H%M%S'`
  $script_dir/pipeline_logger.py \
  -r `cat ${tmp_dir}/nextflow_run_name.txt` \
  -n \${SAMPLE_NAME} \
  -p \${PROCESS_BLOCK} \
  -v 'R --version | head -1' \
  -s \${START_TIME} \
  -e \${STOP_TIME} \
  -d ${log_dir}

  touch stop_flag
  """
}


bcl2fastq_fastqsStopFlagOutChannel
  .concat( fastqcLanesStopFlagOutChannel,
           fastqcSampleStopFlagOutChannel,
           demux_dashStopFlagOutChannel )
  .last()
  .into { logDistillFlagInputChannelCopy01;
          logDistillFlagInputChannelCopy02 }

/*
** Distill the logs for all samples.
*/
process logDistillAllProcess {
  cache 'lenient'
  errorStrategy onError

  input:
    file log_file from logDistillFlagInputChannelCopy01


  script:
  """
  set -ueo pipefail

  ${script_dir}/log_distiller.py -p atac_demux -d ${log_dir} -o ${log_dir}/log.all_samples.txt
  ${script_dir}/log_distiller.py -p atac_demux -d ${log_dir} -s lane pipeline -o ${log_dir}/log.lane.txt
  """
}


Channel
  .fromList( sampleList )
  .set { logDistillSamplesInputChannel }

/*
** Distill individual sample logs.
*/
process logDistillSampleProcess {
  cache 'lenient'
  errorStrategy onError

  input:
    file log_file from logDistillFlagInputChannelCopy01
    val sample_name from logDistillSamplesInputChannel    

  script:
  """
  set -ueo pipefail

  ${script_dir}/log_distiller.py -p atac_demux -d ${log_dir} -s ${sample_name} pipeline -o ${log_dir}/log.${sample_name}.txt
  """
}


/*
** ================================================================================
** Start of Groovy support functions.
** ================================================================================
*/


def printErr( errString ) {
    System.err.println( errString )
}


/*
** Write help information.
*/
def writeHelp() {
    log.info ''
    log.info 'BBI sci-ATAC-seq Demultiplexer'
    log.info '--------------------------------'
    log.info ''
    log.info 'For reproducibility, please specify all parameters in a config file'
    log.info 'and running'
    log.info '  nextflow run main.cf -c CONFIG_FILE.config.'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow run main.nf -c <CONFIG_FILE_NAME>'
    log.info ''
    log.info 'Help: '
    log.info '    --help                              Show this message and exit.'
    log.info ''
    log.info 'Required parameters (specify in your config file):'
    log.info '    params.run_dir = RUN_DIRECTORY             Path to the sequencer run directory.'
    log.info '    params.output_dir = OUTPUT DIRECTORY       Processing output directory.'
    log.info '    params.sample_sheet = SAMPLE_SHEET_PATH    Sample sheet of the format described in the README.'
    log.info ''
    log.info 'Optional parameters (specify in your experiment.config file):'
    log.info '    params.bcl2fastq_cpus = 16                 The number of cores to use for the bcl2fastq run.'
    log.info '    params.max_mem_bcl2fastq = 40              The maximum number of GB of RAM to allocate for bcl2fastq run'
    log.info '    params.index_recipe = 0                    Set explicitly the index recipe. Default is 0 for implicit selection'
    log.info '    params.demux_buffer_blocks = 16            The number of 8K blocks to use for demux output buffer.'
    log.info '    process.maxForks = 20                      The maximum number of processes to run at the same time on the cluster.'
    log.info '    process.queue = "trapnell-short.q"         The queue on the cluster where the jobs should be submitted. '
    log.info ''
    log.info 'Issues? Contact bge@uw.edu'
}


/*
** Report run parameters
*/
def reportRunParams( params, sampleSheetMap ) {

    String s = ""
    s += String.format( "Run parameters\n" )
    s += String.format( "--------------\n" )
    s += String.format( "Sequencing data directory:     %s\n", params.run_dir )
    s += String.format( "Processing output directory:   %s\n", params.output_dir )
    s += String.format( "Processing demux directory:    %s\n", demux_dir )
    s += String.format( "Launch directory:              %s\n", workflow.launchDir )
    s += String.format( "Work directory:                %s\n", workflow.workDir )
    s += String.format( "Sample sheet file:             %s\n", params.sample_sheet )
    s += String.format( "Genomes json file:             %s\n", params.genomes_json )
    s += String.format( "Level:                         %d\n", sampleSheetMap['level'] )
    s += String.format( "Number of wells:               %d\n", sampleSheetMap['number_wells'] )
    s += String.format( "TN5 barcodes:                  %b\n", sampleSheetMap['tn5_barcodes'] )
    s += String.format( "Use all barcodes               %b\n", sampleSheetMap['use_all_barcodes'] )
    s += String.format( "Specify index recipe:          %d\n", params.index_recipe )
    s += String.format( "Maximum bcl2fastq cpus:        %d\n", params.bcl2fastq_cpus )
    s += String.format( "Maximum memory for bcl2fastq:  %d\n", params.max_mem_bcl2fastq )
    s += String.format( "Demux buffer blocks:           %d\n", params.demux_buffer_blocks )
    s += String.format( "\n" )
    print( s )

    /*
    ** Let the operator review and accept or reject the run at this point.
    ** This might be irritating to the operator. If so, comment it out.
    */
    def answer = System.console().readLine 'Do you want to continue ([n]/y)? '
    if( answer != 'y' ) {
        System.exit( 0 )
    }
}


/*
** Check that current version is acceptable given a minimum required version.
*/
def checkVersion( test_version, minimum_version ) {

  def List test_tokens = test_version.tokenize( '.' )
  def List minimum_tokens = minimum_version.tokenize( '.' )

  def num_tokens = Math.min( test_tokens.size(), minimum_tokens.size() )

  def test_ok = true
  for( int i = 0; i < num_tokens; ++i ) {
    def test_int = test_tokens[i].toInteger()
    def minimum_int = minimum_tokens[i].toInteger()
    if( test_int > minimum_int ) {
      break
    }
    if( test_int < minimum_int ) {
      test_ok = false
      break
    }
  }

  return( test_ok )
}


def checkSamplesheetVersion( testVersion, minimumVersion ) {
  if( testVersion == null ) {
    return( false )
  }
  return( checkVersion( testVersion, minimumVersion ) )
}


/*
** Make directory, if it does not exist.
*/
def makeDirectory( directoryName ) {
	dirHandle = new File( directoryName )
	if( !dirHandle.exists() ) {
		if( !dirHandle.mkdirs() ) {
			printErr( "Error: unable to create directory $directoryName" )
			System.exit( -1 )
		}
	}
}


def checkDirectories( params, log_dir, tmp_dir ) {
    /*
    ** Check for existence of run_dir.
    */
    def dhRunDir = new File( params.run_dir )
    if( !dhRunDir.exists() ) {
    	printErr( "Error: unable to find Illumina run directory $params.run_dir" )
    	System.exit( -1 )
    }
    if( !dhRunDir.canRead() ) {
    	printErr( "Error: unable to read Illumina run directory $params.run_dir" )
    	System.exit( -1 )
    }

    /*
    ** Check that either the demux_dir exists or we can create it.
    */
    makeDirectory( demux_dir )

    /*
    ** Check that either the log_dir exists or we can create it.
    */
    makeDirectory( log_dir )

    /*
    ** Check that either the tmp_dir exists or we can create it.
    */
    makeDirectory( tmp_dir )
}


def checkSamplesheet( params ) {
    /*
    ** Check for existence of sample sheet.
    ** Notes:
    **   o  check sample sheet file content
    */
    def fhSampleSheet = new File( params.sample_sheet )
    if( !fhSampleSheet.exists() ) {
    	printErr( "Error: unable to find samplesheet file ${params.sample_sheet}" )
    	System.exit( -1 )
    }
	if( !fhSampleSheet.canRead() ) {
    	printErr( "Error: unable to read file sample sheet ${params.sample_sheet}" )
    	System.exit( -1 )
    }    
}


def archiveRunFiles( params, timeNow )
{
	makeDirectory( "${demux_dir}/reports" )
	file_suffix = timeNow.format( 'yyyy-MM-dd_HH-mm-ss' )
	Path src = Paths.get( params.sample_sheet )
	def ftmp = new File( params.sample_sheet )
	Path dst =  Paths.get( "${demux_dir}/reports/${ftmp.getName()}.${file_suffix}" )
	Files.copy( src, dst )
	def i = 1
	workflow.configFiles.each { aFile ->
		src = aFile
		dst = Paths.get( "${demux_dir}/reports/${aFile.getName()}.${file_suffix}.${i}" )
		Files.copy( src, dst )
		i += 1
	}
}


def readSampleSheetJson( params ) {
    def jsonSlurper = new JsonSlurper()
    sampleSheetMap = jsonSlurper.parse(new File(params.sample_sheet))
    return( sampleSheetMap )
}


def makeSampleList(sampleData) {
  def samples = []

  sampleData['sample_index_list'].each { aSample ->
    samples.add( aSample['sample_id'] )
  }
  samples = samples.unique(mutate=false)
  return(samples)
}


def readIlluminaRunInfo( params ) {
    def command = "${script_dir}/read_run_info.py ${params.run_dir} ATAC-seq"
    def strOut = new StringBuffer()
    def strErr = new StringBuffer()
    def proc = command.execute()
    def jsonSlurper = new JsonSlurper()

    proc.consumeProcessOutput(strOut, strErr)
    proc.waitForProcessOutput()
    if( proc.exitValue() != 0 ) {
        System.err << strErr.toString()
        System.exit( -1 )
    }
    illuminaRunInfoMap = jsonSlurper.parseText(strOut.toString())
    return( illuminaRunInfoMap )
}


def writeArgsJson( params, timeNow ) {
    def i
    
    /*
    ** Set up a Groovy 'Map' called 'mapRunInfo', which has
    **    o  NextFlow runtime parameters
    **    o  sample names and reference genomes
    **    o  Illumina runinfo
    ** Write the information as a json file in demux_dir.
    */
    def mapRunInfo = [:]
    
    /*
    ** Add demux run parameters.
    */
    demuxDict = [:]
    demuxDict['run_date'] = timeNow.format( 'yyyy-MM-dd_HH-mm-ss' )
    demuxDict['run_dir'] = params.run_dir
    demuxDict['output_dir'] = output_dir
    demuxDict['demux_dir'] = demux_dir
    demuxDict['sample_sheet'] = params.sample_sheet
    demuxDict['num_well'] = sampleSheetMap['number_wells']
    demuxDict['level'] = sampleSheetMap['level']
    demuxDict['params.bcl2fastq_cpus'] = params.bcl2fastq_cpus
    demuxDict['max_mem_bcl2fastq'] = params.max_mem_bcl2fastq

    /*
    ** Add sample sheet information.
    */
    def genomeInfo = [:]
    def samples = []
    def peakGroups = [:]
    def peakFiles = [:]
    def fhSampleSheet = new File( params.sample_sheet )
    def jsonSlurper = new JsonSlurper()
    def sampleData = jsonSlurper.parse( fhSampleSheet )

    sampleData['sample_index_list'].each { aSample ->
        samples.add( aSample['sample_id'] )
        genomeInfo.put( aSample['sample_id'], aSample['genome'] )
        peakGroups.put( aSample['sample_id'], aSample['peak_group'] )
        peakFiles.put( aSample['sample_id'], aSample['peak_file'] )
    }


    mapRunInfo['DEMUX'] = demuxDict
    mapRunInfo['sample_data'] = sampleData
    mapRunInfo['samples'] = samples.unique(mutate=false)
    mapRunInfo['genomes'] = genomeInfo
    mapRunInfo['peak_groups'] = peakGroups
    mapRunInfo['peak_files'] = peakFiles

    /*
    ** Add Illumina sequencing run parameters.
    */
    mapRunInfo.put( "sequencing_run_metadata", illuminaRunInfoMap )

    /*
    ** Enclose the mapRunInfo in a map with entry 'RUN001'. This is intended
    ** to simplify aggregating demux runs in an analysis while allowing one
    ** to run an analysis directly from a single demux run output directory.
    */
    def mapRun001 = [:]
    mapRun001 = [ "RUN001": mapRunInfo ]

    /*
    ** Write mapRunInfo to args.json file.
    */
    def args_json = demux_dir + "/args.json"
    File file_json = new File( args_json )
    file_json.write( JsonOutput.prettyPrint( JsonOutput.toJson( mapRun001 ) ) )
}


/*
** A closure for 'publishing' a file to a sample-specific sub-directory.
*/
def qualifyFilename( aFile, aSubDir ) {
    def aSample = aFile.split( '-' )[0]
    def qualifiedFileName = aSample + '/' + aSubDir + '/' + aFile
    return( qualifiedFileName )
}

