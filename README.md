# bbi-sciatac-demux and bbi-sciatac-analyze

*bbi-sciatac-demux* implements Andrew Hill's sci-ATAC-seq demultiplexing processing pipeline and *bbi-sciatac-analyze* implements his sci-ATAC-seq analysis pipeline in the Nextflow domain-specific language.

## Requirements

In order to run these pipeline scripts, Nextflow must be installed on the system where they will run. You can find Nextflow [installation instructions](https://www.nextflow.io/docs/latest/getstarted.html#installation) on the [Nextflow web site](https://www.nextflow.io).

These scripts are written to run on the UW Genome Sciences computing cluster: they depend on the Univa Grid Engine and various installed modules.

Install the scripts using the *git* program to clone the *bbi-lab/bbi-sciatac-demux* and *bbi-lab/bbi-sciatac-analyze* repositories on your computer. Use the commands

```
git clone https://github.com/bbi-lab/bbi-sciatac-demux
git clone https://github.com/bbi-lab/bbi-sciatac-analyze
```

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

The above is a one-time installation setup, or may be required if you need to update the environment.

The bbi-sciatac-analyze pipeline uses two compiled programs to find and process hash reads for sciPlex experiments. These programs require a Rust compiler, which is installed using the information at

https://www.rust-lang.org/tools/install

After installing the Rust compiler, build and install the programs by running the script *bbi-sciatac-analyze/install_rust_programs.sh*.

In addition, a python 3 interpreter is required.


## Overview

#### Nextflow

Nextflow is both a language specification for writing scientific processing pipeline scripts and a program for running such scripts. Important features offered by Nextflow include support for distributed computing on the Univa Grid Engine and the ability to resume processing at a failed processing stage.

The Nextflow pipeline components required for the bbi-sciatac pipelines consist of the

* nextflow program, the executable in the Nextflow distribution,
* nextflow.config file, gives processing parameters, for example, for the Univa Grid Engine,
* params.config file, gives parameters for a particular processing run,
* samplesheet.json file, gives information about the samples in the sequencing run,
* work directory, where Nextflow 'stages' the processing steps. In particular, files created during the processing are stored in sub-directories of work,
* genome file sets,
* computing resources.

#### bbi-sciatac pipeline

The bbi-sciatac pipelines consist of a pair of Nextflow scripts. The first script, *bbi-sciatac-demux/main.nf*, converts Illumina bcl files to fastq files, corrects barcode errors, partitions reads into fastq files by sample, and trims off adapter sequence from the reads. The second script, *bbi-sciatac-analyze/main.nf*, continues processing with read alignments.

#### Samplesheet file

The user-generated samplesheet file (front-end samplesheet) is in CSV format. The program *bbi-sciatac-demux/samplesheet/sciatac_samplesheet.py* converts this file to a JSON file (back-end samplesheet), which is required by the pipeline. In addition, the user can set various processing parameters using *sciatac_samplesheet.py* command line parameters. For additional information about *sciatac_samplesheet.py*, including a detailed description of the front-end samplesheet format, run

```
sciatac_samplesheet.py -d
```

In order to process a sciPlex experiment, run the samplesheet conversion program, *bbi-sciatac-demux/samplesheet/sciatac_samplesheet.py* with the *--hash_file* command line option. The *--hash_file* option takes the full path to the hash index file, and enables hash read processing in the *bbi-sciatac-analyze* pipeline.

#### Genome files

The required genome-related files are specified in the *bbi-sciatac-analyze/genomes.json* file. The file *genomes_format.txt* in the *bbi-sciatac-analyze/genomes* directory describes the *genomes.json* file format. This description is from Andrew Hill.

#### nextflow.config file

The *nextflow.config* file defines processing values such as the required modules, memory, and number of CPUs for each processing stage, which do not change typically from run-to-run. The file can be left in the bbi-sciatac-\* installation directory where Nextflow searches for it automatically when a pipeline run starts. The supplied *nextflow.config* file has two profiles: the default profile, called *standard*, defines modules used by the pipeline on CentOS 7 systems in the UW Genome Sciences cluster, and the *centos_6* profile, which defines modules used by the pipeline on CentOS 6 systems in the UW Genome Sciences cluster. In order to run the pipelines with the *centos_6* profile, add the command line parameter `-profile centos_6` to the nextflow run command, for example


```
nextflow run bbi-dmux -profile centos_6 -c experiment.config
```

This *nextflow.config* file has comments that give additional information.

#### Parameters configuration file

The lines in the *experiment.config* file define parameters required by the processing pipeline, such as the Illumina run directory and the processing pipeline output directory. There is additional information in the example *experiment.config* files in the pipeline repositories.

#### Running the pipeline scripts.

Run the pipeline scripts in a cluster node a significant memory resource. For example, use

```
qlogin -l mfree=16G
```

Copy the *experiment.config* file and bash scripts *bbi-sciatac-demux/scripts/run.sciatac-demux.sh* and *bbi-sciatac-analyze/scripts/run.sciatac-analyze.sh* to the directory where you want the pipeline output files. Edit *run.sciatac-demux.sh* and *n.sciatac-analyze.sh*as described in the script comments. Edit *experiment.config* as required then run *run.sciatac-demux.sh* followed by *run.sciatac-analyze.sh*.

```
bash run.sciatac-demux.sh
.
.
.
bash run.sciatac-analyze.sh
```

Each script starts by reporting important parameters and then waits run confirmation.

The script *run.sciatac-demux.sh* creates a sub-directory called *demux_out* where the demux pipeline writes its output files. The *run.sciatac-analyze.sh* script creates a sub-directory called *analyze_out* where the analyze pipeline writes its output files.

#### Additional information

We advise users to run the pipelines in tmux sessions so that closing the terminal window does not terminate the pipeline.

When run with the *-with-trace <trace_file_name>* command line parameter, Nextflow writes helpful information to a trace file. In the case that Nextflow is run on a Grid Engine, the trace file has a column labelled *native_id*, which is the Grid Engine job id. You can use this job id with the *qacct -j <job_id>* Grid Engine command to better understand why the Grid Engine may have killed a process. Additionally, the trace file has columns with the amount of memory used by each process and process hash strings, which you can use to tune the Grid Engine memory requests and find the corresponding process sub-directory in the Nextflow work directory, respectively.

#### Questions and errors:
If you run into problems, please leave a detailed description in the issue tab above.

## Acknowledgements

A debt of gratitude is due to the many members of the Shendure and Trapnell labs as well as the BBI team who have worked on portions of this pipeline or its predecessors, especially Hannah Pliner, Aishwarya Gogate, Andrew Hill, and Jonathan Packer.
