#!/bin/bash

module load bcl2fastq/2.20
module load kentUtils/302.1.0
module load tabix/0.2.6

script="./pipeline_logger.py"

run_name='fully_here'
sample_name='NA'
process_block='01_test_block'
start_time='20210420:160431'
stop_time='20210420:160445'
log_dir='log_dir'

value1='value1'
value2='value2'

comm1="test_comm $value1"
comm2="test_comm $value1 $value2"

${script} \
-r ${run_name} \
-n ${sample_name} \
-p ${process_block} \
-v 'awk --version | head -1' 'zcat --version | head -1' 'bedToBigBed 2>&1 > /dev/null | head -1' 'bgzip --version | head -1' 'tabix 2>&1 > /dev/null | head -3 | tail -1' 'zcat --version | head -1' \
-s ${start_time} \
-e ${stop_time} \
-f "ATAC.Sentinel-RUN001_L001.20210504143840.adapter_trimming.47_f11c6f77fd1aab74e56c166d3cf615.log" \
-d ${log_dir} \
-c "$comm1" "$comm2"
