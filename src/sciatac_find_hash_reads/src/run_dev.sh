#!/bin/bash

DATA_DIR="/home/brent/git/bbi-hash-experiment-scripts/src/atac-seq/bge/data_test"


echo "Run on uncompressed 1M read fastq files."
cargo run -- -1 $DATA_DIR/r1.fq -2 $DATA_DIR/r2.fq -h $DATA_DIR/hash.tsv -o hits.out

