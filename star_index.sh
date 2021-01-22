#! /bin/bash

### This script creates the index for the STAR alignment method.

[[ $# -ne 2 ]] && echo "Usage: ./star_index.sh genome.fasta nthreads" && exit 1

fasta=$1
n=$2
build=$(basename $fasta | cut -f 1 -d '.')
index_dir=/project/MDL_ChIPseq/process_sequencing_data/genome/star_index
fasta_dir=/project/MDL_ChIPseq/process_sequencing_data/genome/fasta

[[ ! -e $fasta_dir/$build.fa ]] && ln -s $fasta $fasta_dir/$build.fa
fasta=$fasta_dir/$build.fa 

STAR --runThreadN $n \
	 --runMode genomeGenerate \
	 --genomeDir $index_dir/${build} \
	 --genomeFastaFiles $fasta
