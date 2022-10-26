#!/bin/bash

#
# qlogin -l mfree=8G
#

#
# Merge trimmed fastq files of two (or more) bbi-sciatac-demux runs into fastq_trim directories.
# Set up:
#   o  list names of samples for which to merge fastq files in sample_list
#   o  set root directory of first source demux directory in src1_root_dir
#   o  set root directory of second source demux directory in src2_root_dir
#   o  set root of merged demux directory in dst_root_dir
#
# The log file has line counts of the input and output files.
#
# Check the output with
#   cat run_merge.log4 | grep '^Line' | awk '{if($7!=$9){print$7,$8,$9}}'
#

log_file='run_merge.log'


#sample_list="W134.aorta W134.inferior.vena.cava W134.spleen W136.lung W136.spleen W137.liver W137.lung W137.pancreas W137.spleen W139.lung W139.pancreas W139.spleen W141.lung W141.pancreas W141.spleen W142.aorta W142.inferior.vena.cava W142.lung W142.spleen W158.LV.heart W135.pancreas NSM198 CH12.LX Barn"

sample_list="NSM198"

src1_root_dir="/net/bbi/vol1/data/sciATACseq/novaseq_runs/ATAC3-004-a"
src2_root_dir="/net/bbi/vol1/data/sciATACseq/novaseq_runs/ATAC3-005-a"
dst_root_dir="/net/bbi/vol1/data/sciATACseq/novaseq_runs/ATAC3-004_005-a.r3"

echo "Check that untrimmed and trimmed fastq files exist in both source directories."
errorFlag=0
for sample in $sample_list
do
  echo $sample

  fastq_dir1="${src1_root_dir}/demux_out/${sample}/fastqs_barcode" 
  fastq_dir2="${src2_root_dir}/demux_out/${sample}/fastqs_barcode" 
  lfil1=`pushd ${fastq_dir1} > /dev/null; ls *.gz; popd > /dev/null` 
  lfil2=`pushd ${fastq_dir2} > /dev/null; ls *.gz; popd > /dev/null` 
  for fil in $lfil1 
  do 
     echo $fil 
  done 

  if [ "$lfil1" != "$lfil2" ] 
  then 
    echo "Error: fastqs_barcode file lists not equal for sample ${sample}" 
    errorFlag=1 
  fi 


  fastq_dir1="${src1_root_dir}/demux_out/${sample}/fastqs_trim"
  fastq_dir2="${src2_root_dir}/demux_out/${sample}/fastqs_trim"
  lfil1=`pushd ${fastq_dir1} > /dev/null; ls *.gz; popd > /dev/null`
  lfil2=`pushd ${fastq_dir2} > /dev/null; ls *.gz; popd > /dev/null`
  for fil in $lfil1
  do
     echo $fil
  done

  if [ "$lfil1" != "$lfil2" ]
  then
    echo "Error: file fastqs_trim lists not equal for sample ${sample}"
    errorFlag=1
  fi
done

if [ "$errorFlag" -eq 1 ]
then
  exit -1
fi

echo "Make merged output directory ${dst_root_dir}."
mkdir -p ${dst_root_dir}/demux_out

echo "Make merged output sample fastq_barcode directories." 
for sample in $sample_list 
do 
  echo $sample 
  fastq_dir_merged="${dst_root_dir}/demux_out/${sample}/fastqs_barcode" 
  mkdir -p ${fastq_dir_merged} 
done 

#echo "Make merged output sample fastq_trim directories."
#for sample in $sample_list
#do
#  echo $sample
#  fastq_dir_merged="${dst_root_dir}/demux_out/${sample}/fastqs_trim"
#  mkdir -p ${fastq_dir_merged}
#done


echo "Merge barcode fastq files." 
for sample in $sample_list 
do 
  echo "Process sample $sample" 
 
  fastq_dir1="${src1_root_dir}/demux_out/${sample}/fastqs_barcode" 
  fastq_dir2="${src2_root_dir}/demux_out/${sample}/fastqs_barcode" 
  lfil=`pushd ${fastq_dir1} > /dev/null; ls *.gz; popd > /dev/null` 
  fastq_dir_merged="${dst_root_dir}/demux_out/${sample}/fastqs_barcode" 
 
  pushd ${fastq_dir_merged} 
  for fil in $lfil 
  do 
    echo "Merge files $fil" 
 
    cat ${fastq_dir1}/${fil} ${fastq_dir2}/${fil} > ${fastq_dir_merged}/$fil

    echo "Check files $fil" 

    echo "-----" >> ${dst_root_dir}/${log_file} 
    c1=`zcat ${fastq_dir1}/${fil} | wc -l` 
    c2=`zcat ${fastq_dir2}/${fil} | wc -l` 
    c3=`zcat ${fastq_dir_merged}/${fil} | wc -l` 
    echo "File: $fil" >> ${dst_root_dir}/${log_file} 
    echo "Line counts: $c1 + $c2 = $((c1+c2)) ?= $c3" >> ${dst_root_dir}/${log_file} 
  done 
  popd 

done 


#echo "Merge trimmed fastq files."
#for sample in $sample_list
#do
#  echo "Process sample $sample"
#
#  fastq_dir1="${src1_root_dir}/demux_out/${sample}/fastqs_trim"
#  fastq_dir2="${src2_root_dir}/demux_out/${sample}/fastqs_trim"
#  lfil=`pushd ${fastq_dir1} > /dev/null; ls *.gz; popd > /dev/null`
#  fastq_dir_merged="${dst_root_dir}/demux_out/${sample}/fastqs_trim"
#
#  pushd ${fastq_dir_merged}
#  for fil in $lfil
#  do
#    echo "Merge files $fil"
#
#    cat ${fastq_dir1}/${fil} ${fastq_dir2}/${fil} > ${fastq_dir_merged}/$fil
#
#    echo "Check files $fil"
#
#    echo "-----" >> ${dst_root_dir}/${log_file}
#    c1=`zcat ${fastq_dir1}/${fil} | wc -l`
#    c2=`zcat ${fastq_dir2}/${fil} | wc -l`
#    c3=`zcat ${fastq_dir_merged}/${fil} | wc -l`
#    echo "File: $fil" >> ${dst_root_dir}/${log_file}
#    echo "Line counts: $c1 + $c2 = $((c1+c2)) ?= $c3" >> ${dst_root_dir}/${log_file}
#  done
#  popd
#
#done

