# bbi-sciatac-demux and bbi-sciatac-analyze

*bbi-sciatac-demux* implements Andrew Hill's sci-ATAC-seq demultiplexing processing pipeline in the NextFlow domain-specific language and *bbi-sciatac-analyze* implements Andrew Hill's sci-ATAC-seq analysis processing pipeline.

## Requirements

In order to run these scripts, Nextflow must be installed on the system where you will run the processing pipeline. You can find Nextflow [installation instructions](https://www.nextflow.io/docs/latest/getstarted.html#installation) on the [Nextflow web site](https://www.nextflow.io).

These scripts are written to run on the UW Genome Sciences computing cluster: they depend on the Univa Grid Engine and various installed modules.

For efficiency and convenience, they are written to run certain stages in python virtual environments. The environments may need to be built on a node with the CPU architecture on which you will run the scripts. In order to build the environments use the commands

```
qlogin -l mfree=16G
module load git/latest

cd bbi-sciatac-demux
bash create_virtual_envs.sh
cd ..

cd bbi-sciatac-analyze
bash create_virtual_envs.sh
cd ..
```

The above is a one-time installation setup, or may be required if you need to update the pipeline.

## Overview

#### Nextflow

Nextflow is both a language specification for writing scientific processing pipeline scripts and a program for running such scripts. Important features offered by Nextflow include support for distributed computing on the Univa Grid Engine and the ability to resume processing at a failed stage processing step.

The Nextflow pipeline components required for the bbi-sciatac pipeline consist of the

* nextflow program, the executable in the Nextflow distribution,
* nextflow.config file, gives processing parameters, for example, the Univa Grid Engine,
* params.config file, gives parameters for a particular processing run,
* work directory, where Nextflow 'stages' the processing steps. In particular, files created during the processing are stored in sub-directories of work,
* genome file sets,
* computing resources.

#### bbi-sciatac pipeline

The bbi-sciatac pipeline consists of a pair of Nextflow scripts. The first converts an Illumina bcl file to fastq files, corrects barcode errors, partitions reads into fastq files by sample, and trims off adapter sequence from the reads. The second script continues processing with read alignments through to making count matrices.

Note: I no longer support the *setup_sciatac.py* script. Instead, consider using the scripts *bbi-sciatac-demux/run.sciatac-demux.sh* and *bbi-sciatac-analyze/run.sciatac-analyze.sh* I leave *setup_sciatac.py* here because it may have some value as a guide for setting up scripts.

The bbi-sciatac-demux repository includes a script called *setup_sciatac.py*, which sets up the required directories and files for the pipeline runs. The script prompts for required values and allows editing of values. The editable values consist of

* Processing (parent) directory
* Illumina run directory
* Samplesheet file
* Levels
* Number of wells
* Demux script name
* Demux output directory
* Analyze script name
* Analyze output directory
* Nextflow parameters filename (for example, params.config)
* Maximum number of CPUs
* Maximum memory for bcl2fastq
* Grid engine queue name (optional)
* Grid engine cluster options (optional)
* Genomes JSON file

Additionally, you will need to set the following values in the script

* nextflowExe, the path to the Nextflow executable,
* demuxNextflowScript, the path to the demux Nextflow script, which is called main.nf,
* analyzeNextflowScript, the path to the analyze Nextflow script, which is called main.nf.

#### Samplesheet file

The samplesheet is a tab-delimited file with the format:

```
sample_id   ranges  genome
sample_1    1-96:25-48:65-80:1-20   hg19
sample_2    1-96:25-48:65-80:21-40  hg19
sample_3    1-96:25-48:65-80:41-60  mm9
sample_4    1-96:25-48:65-80:61-80  hg19_mm9
sample_5    1-96:25-48:65-80:80-96  mm9
```

The `sample_id` column can have whatever you want as an ID for each sample in the run.

The `ranges` column specifies the range of `<N7 ligation indices>:<N7 PCR indices>:<N5 PCR indices>:<N5 ligation indices>` used for this run. Note that the N7/P7 indices are numbered by column and the N5/N7 indices are numbered by row. Multiple ranges for a barcode are separated by commas.

There is a rudimentary script for converting a BBI CSV samplesheet to the required format in the *bbi-sciatac-demux/samplesheet* directory.

Warning: sample names cannot have dashes in them, or any special character that does not appear typically in a Unix file name so such characters are converted to dots (any character other than 'a-z', 'A-Z', '0-9', '_',  and '.'). This is because the script adds the sample name to the main file name, separated by a dash, so that it can identify the sample based the filename. For example, the file 'sample1-RUN001_L001_R1.fastq.gz' is from a sample called 'sample1'.

#### Genome files

The required genome-related files are specified in the *bbi-sciatac-analyze/genomes.json* file. The required files include

* bowtie genome indexes
* genome whitelist regions
* chromosome sizes

The *bbi-sciatac-analyze* repository has some scripts for assisting with building these files in the *genomes* subdirectory. These scripts are from Andrew Hill.

#### nextflow.config file

#### Parameters configuration file

#### Demultiplexing parameters

#### Analyze parameters

#### Running the pipeline scripts.

Run the pipeline scripts in a cluster node with significant a significant memory resource. For example, use

```
qlogin -l mfree=16G
```

Run the bash script *run.demux.sh* first. After that finishes, run *run.analyze.sh*.

```
bash run.demux.sh
.
.
.
bash run.analyze.sh
```

The *bbi-sciatac-demux* pipeline creates the file *args.json* in the demux output directory. The *bbi-sciatac-analyze* pipeline script looks for *args.json* in the demux output directory, and uses it to read the trimmed fastq files from the *<sample_name>/fastqs_trim* sub-directory in the demux output directory, and continues the analysis starting with read alignments.

#### Additional information

We advise users to run the pipelines in tmux sessions so closing the terminal window does not terminate the pipeline.

#### Questions and errors:
If you run into problems, please leave a detailed description in the issue tab above.

## Acknowledgements

A debt of gratitude is due to the many members of the Shendure and Trapnell labs as well as the BBI team who have worked on portions of this pipeline or its predecessors, especially Hannah Pliner, Aishwarya Gogate, Andrew Hill, and Jonathan Packer.
