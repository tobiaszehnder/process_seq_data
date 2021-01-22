#! /bin/bash

[[ $# < 4 ]] && echo "Usage: ./star_align.sh build nthreads outfile fastq[both for paired-end]" && exit 1

build=$1
n=$2
outfile=$3 # path/to/file.bam
outfile_prefix=${outfile%.*} # /path/to/file
outdir=${outfile%/*} # /path/to
sample_name=${outfile_prefix##*/} # file
log_dir=${outdir}/log
mkdir -p $log_dir
nfastq=$(($#-3))
fastq="${@:4:nfastq}"
star_index_dir=/project/MDL_ChIPseq/process_sequencing_data/genome/star_index/$build
star_index_dir=/project/bat_analysis/Software/STAR/Indices/carPer2
tmp=/scratch/local2/${sample_name}_tmp
[[ -e $tmp ]] && rm -rf $tmp # clear tmp folder if it exists

# STAR=STAR
STAR=/scratch/ngsvin/bin/mapping/STAR/STAR_2.6.1d/bin/Linux_x86_64_static/STAR

$STAR --runMode alignReads \
	  --genomeLoad LoadAndKeep \
	  --genomeDir $star_index_dir \
	  --runThreadN $n \
	  --readFilesIn $fastq \
	  --alignEndsType EndToEnd \
	  --alignIntronMax 1 \
	  --alignSJDBoverhangMin 999 \
	  --limitBAMsortRAM 20000000000 \
	  --outFileNamePrefix ${log_dir}/${sample_name}. \
	  --outSAMtype BAM SortedByCoordinate \
	  --outSAMattributes Standard \
	  --outSAMunmapped None \
	  --outSAMmode NoQS \
	  --outSAMattrRGline ID:$sample_name \
	  --outStd BAM_SortedByCoordinate \
	  --outFilterMultimapNmax 5 \
	  --outFilterMatchNminOverLread 0.94 \
	  --outTmpDir $tmp \
	  --readFilesCommand 'gzip -dc' \
	  > $outfile
