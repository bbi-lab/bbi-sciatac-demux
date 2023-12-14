#!/bin/bash

PGRM="/home/brent/git/bbi-hash-experiment-scripts/src/atac-seq/bge/sciatac_find_hash_reads/target/release/sciatac_find_hash_reads"

DATA_DIR="/home/brent/git/bbi-hash-experiment-scripts/src/atac-seq/bge/sciatac_find_hash_reads/data_test"


echo "Run on uncompressed 10M read fastq files."
time $PGRM -1 $DATA_DIR/r1.10m.fq -2 $DATA_DIR/r2.10m.fq -h $DATA_DIR/hash.tsv -o hits.out > /dev/null

echo "Run on compressed 10M read fastq files."
time $PGRM -1 $DATA_DIR/r1.10m.fq.gz -2 $DATA_DIR/r2.10m.fq.gz -h $DATA_DIR/hash.tsv -o hits.out > /dev/null


