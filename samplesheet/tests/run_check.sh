#!/bin/bash


lfil="ATAC3-018_samplesheet samplesheet.test01 samplesheet.test02 samplesheet.test03 samplesheet.test04 samplesheet.test05 samplesheet.test06 samplesheet.test07 samplesheet.test08 samplesheet.test09 samplesheet.test10 sciatac_samplesheet_example_01 sciatac_samplesheet_example_02 sciatac_samplesheet_example_07 sciatac_samplesheet_example_08"

#lfil="ATAC3-018_samplesheet"

PGRM1="../sciatac_samplesheet.py"
PGRM2="../check_sciatac_samplesheet.py"
OUT_DIR="out_files.20231116"

for fil in $lfil
do
  echo "${fil}.csv"
  ${PGRM1} -i ${fil}.csv -o "${OUT_DIR}/${fil}.json"
  ${PGRM2} -i "${OUT_DIR}/${fil}.json" > "${OUT_DIR}/${fil}.check"
done
