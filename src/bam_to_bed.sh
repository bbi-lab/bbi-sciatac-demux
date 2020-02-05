#!/usr/bin/env bash
module purge
module load modules modules-init modules-gs
module load bedtools/2.26.0

BAM=$1;
BED=$2;

bedtools bamtobed -i $BAM | cut -f1-3 | gzip > $BED
