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

import groovy.json.JsonOutput
import java.nio.file.Files
import java.nio.file.Path 
import java.nio.file.Paths

/*
** Run date/time.
*/
def timeNow = new Date()

/*
** Where to find scripts.
** Note: script_dir needs to be visible within Groovy functions
**       so there is no 'def', which makes it global.
**
*/
def pipeline_path="/net/gs/vol1/home/bge/git/bbi-sciatac-demux"
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
*/
params.help = false
params.num_well = 384
params.level = 3
params.max_cores = 16
params.max_mem_bcl2fastq = 40

/*
** Initialize optional parameters to null.
*/
// none at this time

/*
** Define/initialize some internal parameters.
** Notes:
**   o  sample_sheet_phantom is a non-existent file required by bcl2fastq
*/
def sample_sheet = params.sample_sheet
def run_dir      = params.run_dir
def demux_dir    = params.demux_dir
def num_threads_bcl2fasta_process
def mem_bcl2fastq = params.max_mem_bcl2fastq
def num_threads_bcl2fastq_io
def options_barcode_correct = ''
def sample_sheet_phantom = run_dir + '/Sample_sheet.csv'

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
assert ( params.run_dir && params.demux_dir && params.sample_sheet && params.num_well && params.level ) : "missing config file: use -c CONFIG_FILE.config that includes run_dir, demux_dir, sample_sheet, num_well, and level parameters"
assert ( params.num_well == 96 || params.num_well == 384 ) : "invalid params.num_well value: must be either 96 or 384"
assert ( params.level == 2 || params.level == 3 ) : "invalid params.level value: must be either 2 or 3"
assert ( ( params.level == 2 && params.num_well == 96 ) || ( params.level == 3 && params.num_well == 96 ) || ( params.level == 3 && params.num_well == 384 ) ) : "invalid combination of params.level and params.num_well"

/*
** Report run parameter values.
*/
reportRunParams( params )

/*
** Archive configuration and samplesheet files in demux_dir.
*/
archiveRunFiles( params, timeNow )

/*
** Check that required directories exist or can be made.
*/
checkDirectories( params )

/*
** Check sample sheet file.
*/
checkSamplesheet( params )

/*
** Read Illumina run information.
*/
illuminaRunInfoMap = readIlluminaRunInfo( params )

/*
** Write run information to args.json file.
*/
writeArgsJson( params, timeNow )

/*
** Does the i5 sequence require reverse complementing?
** If yes, set sequence_flag.
*/
def sequencer_flag = ( illuminaRunInfoMap['reverse_complement_i5'].toInteger() ) ? '-X' : ''

/*
** Add optional barcode correction parameters.
*/
if( params.level == 2 && params.num_well == 96 ) {
	options_barcode_correct += ' --two_level_indexed_tn5'
}
if( params.num_well == 384 ) {
	options_barcode_correct += ' --wells_384'
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
if( params.max_cores > 16 ) {
	num_threads_bcl2fasta_process = 16
	mem_bcl2fastq = params.max_mem_bcl2fastq / 16
} else {
	num_threads_bcl2fasta_process = params.max_cores
	mem_bcl2fastq = params.max_mem_bcl2fastq / num_threads_bcl2fasta_process
}

if( num_threads_bcl2fasta_process / 2 < 4 ) {
	num_threads_bcl2fastq_io = num_threads_bcl2fasta_process / 2
} else {
	num_threads_bcl2fastq_io = 4
}

/*
** Run bcl2fastq to make fastq files from Illumina bcl files.
**
** Notes:
**   o  copy 'Reports' and 'Stats' directories to an accessible directory.
**   o  bcl2fastq/2.16 appears to be a special version. (bcl2fastq/2.20 does not write the required barcodes to sequence headers)
**      I understand that the more recent bcl2fastq programs required a sample sheet with Ns as sequences.
*/
Channel
    .fromList( [ "${sample_sheet_phantom}" ] )
    .set { bcl2fastqInChannel }
    
process bcl2fastq {
  cache 'lenient'
  errorStrategy onError
  penv 'serial'
  cpus num_threads_bcl2fasta_process
  memory "$mem_bcl2fastq GB"    
  module 'java/latest:modules:modules-init:modules-gs'
  module 'mpfr/3.1.0:mpc/0.8.2:gmp/5.0.2:gcc/4.9.1:bcl2fastq/2.16'
//  note: bge: the following publishDir statement is commented out in bbi_dmux. I want it to work for diagnostics during development.
  publishDir path: "$params.demux_dir/fastqs_bcl2fastq", pattern: "Undetermined_S0_*.fastq.gz", mode: 'copy'
  publishDir path: "$params.demux_dir/fastqs_bcl2fastq", pattern: "Stats/*", mode: 'copy'
  publishDir path: "$params.demux_dir/fastqs_bcl2fastq", pattern: "Reports/*", mode: 'copy'

  input:
    val inPhantom from bcl2fastqInChannel

  output:
    set file("Undetermined_S0_*_R1_001.fastq.gz"), file("Undetermined_S0_*_R2_001.fastq.gz") into bcl2fastq_fastqs mode flatten
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
    bcl2fastq --runfolder-dir      ${params.run_dir} \
              --output-dir         . \
              --interop-dir        . \
              --sample-sheet       ${sample_sheet_phantom} \
              --loading-threads    ${num_threads_bcl2fastq_io} \
              --processing-threads ${num_threads_bcl2fasta_process}  \
              --writing-threads    ${num_threads_bcl2fastq_io} \
              --ignore-missing-positions \
              --ignore-missing-controls \
              --ignore-missing-filter \
              --ignore-missing-bcls
    """
    //               --tiles s_1 \
}

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
**   o  the barcode_stats_json and barcode_counts_csv output channels are dummy
**      channels. It appears that the files do not get copied by the publishDir
**      directive if these files are not in a channel.
*/
process barcode_correct {
  cache 'lenient'
  errorStrategy onError
  memory "16 GB"
  module 'java/latest:modules:modules-init:modules-gs:zlib/1.2.6:pigz/latest'
  publishDir    path: "${params.demux_dir}", saveAs: { qualifyFilename( it, "fastqs_barcode" ) }, pattern: "*.fastq.gz", mode: 'copy'
  publishDir    path: "${params.demux_dir}/fastqs_barcode", pattern: "*.stats.json", mode: 'copy'
  publishDir    path: "${params.demux_dir}/fastqs_barcode", pattern: "*.barcode_counts.csv", mode: 'copy'
  publishDir    path: "${params.demux_dir}/fastqs_barcode", pattern: "*.index_counts.csv", mode: 'copy'
  publishDir    path: "${params.demux_dir}/fastqs_barcode", pattern: "*.tag_pair_counts.csv", mode: 'copy'
  publishDir    path: "${params.demux_dir}/fastqs_barcode", pattern: "*.pcr_pair_counts.csv", mode: 'copy'
  
  input:
    // set file(R1), file(R2) from bcl2fastq_fastqs.splitFastq(by: params.fastq_chunk_size, file: true, pe: true)
    set file(R1), file(R2) from bcl2fastq_fastqs
    
  output:
    file "*.fastq.gz" into barcode_fastqs mode flatten
    file "*.stats.json" into barcode_stats_json mode flatten
    file "*.barcode_counts.csv" into barcode_counts_csv mode flatten
    file "*.index_counts.csv" into index_counts_csv mode flatten
    file "*.tag_pair_counts.csv" into tag_pair_counts_csv mode flatten
    file "*.pcr_pair_counts.csv" into pcr_pair_counts_csv mode flatten

  script:
  """
  source $pipeline_path/load_pypy_env_reqs.sh
  PS1=\${PS1:-}
  source $script_dir/pypy_env/bin/activate

  pypy $script_dir/barcode_correct_sciatac.py \
                       --samplesheet $sample_sheet \
                       -1 <(zcat $R1) \
                       -2 <(zcat $R2) \
                       --filename $R1 \
                       --out_dir . \
                       --stats_out 1 \
                       $options_barcode_correct \
                       $sequencer_flag
  deactivate
  """
}

/*
** Gather R1/R2 pairs of fastq filenames because
** Trimmomatic processes read pairs using separate
** input/output files for reads 1 and 2.
*/
get_prefix_barcode_fastqs = { file -> (file - ~/_R[12]\.fastq\.gz/) }

barcode_fastqs
  .map { file -> tuple( get_prefix_barcode_fastqs(file.name), file) }
  .groupTuple()
  .set { barcode_fastqs_paired }

/*
** Run adapter trimming script.
*/
def num_threads_trimmomatic = 4
def mem_trimmomatic = 4.0 / num_threads_trimmomatic
def trimmomatic_exe="${script_dir}/Trimmomatic-0.36/trimmomatic-0.36.jar"
def adapters_path="${script_dir}/Trimmomatic-0.36/adapters/NexteraPE-PE.fa:2:30:10:1:true"
process adapter_trimming {
  cache 'lenient'
  errorStrategy onError
  penv 'serial'
  cpus num_threads_trimmomatic
  memory "${mem_trimmomatic} GB"
  module 'java/8u25:modules:modules-init:modules-gs'
  publishDir path: "$params.demux_dir", saveAs: { qualifyFilename( it, "fastqs_trim" ) }, pattern: "*.fastq.gz", mode: 'copy'
  
  input:
  set prefix, file(read_pair) from barcode_fastqs_paired
  
  output:
  set file("*_R1.trimmed.fastq.gz"), file("*_R2.trimmed.fastq.gz"), file("*_R1.trimmed_unpaired.fastq.gz"), file("*_R2.trimmed_unpaired.fastq.gz") into fastqs_trim
  
  script:
  """
  mkdir -p fastqs_trim
  java -Xmx1G -jar $trimmomatic_exe \
       PE \
       -threads $num_threads_trimmomatic \
       $read_pair \
       ${prefix}_R1.trimmed.fastq.gz \
       ${prefix}_R1.trimmed_unpaired.fastq.gz \
       ${prefix}_R2.trimmed.fastq.gz \
       ${prefix}_R2.trimmed_unpaired.fastq.gz \
       ILLUMINACLIP:${adapters_path} \
       TRAILING:3 \
       SLIDINGWINDOW:4:10 \
       MINLEN:20
  """
}


/*
** This is the end of the 'fastq-producing' sciatac_pipeline processing steps.
** However, this script requires still statistics/report generating functions.
*/



/*
** ================================================================================
** Start of Groovy support functions.
** ================================================================================
*/


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
    log.info '    params.demux_dir = OUTPUT DIRECTORY        Processing output directory.'
    log.info '    params.sample_sheet = SAMPLE_SHEET_PATH    Sample sheet of the format described in the README.'
    log.info '    params.num_well = 384                      Number of u-titer plate wells.'
    log.info '    params.level = 3                           Level of run (2 or 3).'
    log.info ''
    log.info 'Optional parameters (specify in your config file):'
    log.info '    params.max_cores = 16                      The maximum number of cores to use - fewer will be used if appropriate.'
    log.info '    process.maxForks = 20                      The maximum number of processes to run at the same time on the cluster.'
    log.info '    process.queue = "trapnell-short.q"         The queue on the cluster where the jobs should be submitted. '
    log.info '    params.max_mem_bcl2fastq = 40              The maximum number of GB of RAM to allocate for bcl2fastq run'
    log.info ''
    log.info 'Notes:'
    log.info '  o  valid combinations of (params.level, params.num_well) are (3, 384), (3, 96), and'
    log.info '     (2, 96). The last combination is for the original barcoded Tn5-based assay.'
    log.info ''
    log.info 'Issues? Contact bge@uw.edu'
}


/*
** Report run parameters
*/
def reportRunParams( params ) {

    String s = ""
    s += String.format( "Run parameters\n" )
    s += String.format( "--------------\n" )
    s += String.format( "Sequencing data directory:     %s\n", params.run_dir )
    s += String.format( "Processing output directory:   %s\n", params.demux_dir )
    s += String.format( "Launch directory:              %s\n", workflow.launchDir )
    s += String.format( "Work directory:                %s\n", workflow.workDir )
    s += String.format( "Combinatorial levels:          %d\n", params.level )
    s += String.format( "Sample sheet file:             %s\n", params.sample_sheet )
    s += String.format( "Maximum cores:                 %d\n", params.max_cores )
    s += String.format( "Maximum memory for bcl2fastq:  %d\n", params.max_mem_bcl2fastq )
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


def checkDirectories( params ) {
    /*
    ** Check for existence of run_dir.
    */
    def dhRunDir = new File( params.run_dir )
    assert dhRunDir.exists() : "unable to find Illumina run directory $run_dir"
    assert dhRunDir.canRead() : "unable to read Illumina run directory $run_dir"

    /*
    ** Check that either the demux_dir exists or we can create it.
    */
    def dhOutDir = new File( params.demux_dir )
    if( !dhOutDir.exists() ) {
       assert dhOutDir.mkdirs() : "unable to create output directory $demux_dir"
    }
}


def checkSamplesheet( params ) {
    /*
    ** Check for existence of sample sheet.
    ** Notes:
    **   o  check sample sheet file content
    */
    def fhSampleSheet = new File( params.sample_sheet )
    assert fhSampleSheet.exists() : "unable to find sample sheet file ${params.sample_sheet}"
    assert fhSampleSheet.canRead() : "unable to read file sample sheet ${params.sample_sheet}"
}


def archiveRunFiles( params, timeNow )
{
  file_suffix = timeNow.format( 'yyyy-MM-dd_HH-mm-ss' )
  Path src = Paths.get( params.sample_sheet )
  def ftmp = new File( params.sample_sheet )
  Path dst =  Paths.get( "${params.demux_dir}/${ftmp.getName()}.${file_suffix}" )
  Files.copy( src, dst )
  def i = 1
  workflow.configFiles.each { aFile ->
    src = aFile
    dst = Paths.get( "${aFile.getName()}.${file_suffix}.${i}" )
    Files.copy( src, dst )
    i += 1
  }
}


def readIlluminaRunInfo( params ) {
    def command = "${script_dir}/run_info_read.py ${params.run_dir}"
    def strOut = new StringBuffer()
    def strErr = new StringBuffer()
    def proc = command.execute()
    proc.consumeProcessOutput(strOut, strErr)
    proc.waitForProcessOutput()
    if( proc.exitValue() != 0 ) {
        System.err << strErr.toString()
        System.exit( -1 )
    }

    def illuminaRunInfoMap = [:]
    def tokens = strOut.tokenize()
    for( i = 0; i < tokens.size; i += 2 ) {
        illuminaRunInfoMap.put( tokens[i], tokens[i+1] )
    }
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
    demuxDict['demux_dir'] = params.demux_dir
    demuxDict['sample_sheet'] = params.sample_sheet
    demuxDict['num_well'] = params.num_well
    demuxDict['level'] = params.level
    demuxDict['params.max_cores'] = params.max_cores
    demuxDict['max_mem_bcl2fastq'] = params.max_mem_bcl2fastq

    /*
    ** Add sample sheet information.
    */
    def genomeInfo = [:]
    def samples = []
    def fhSampleSheet = new File( params.sample_sheet )
    fhSampleSheet.eachLine { line ->
        def (sample, wells, genome) = line.split(/[ \t]+/)
        if( sample != "sample_id" ) {
            samples.add( sample )
            genomeInfo.put( sample, genome )
        }
     }
     mapRunInfo['DEMUX'] = demuxDict
     mapRunInfo['samples'] = samples
     mapRunInfo['genomes'] = genomeInfo

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
    def args_json = params.demux_dir + "/args.json"
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

