#! /bin/bash

set -e # exit when any command fails

usage()
{
cat << EOF
usage: ${0##*/} [options...] [R1.fastq.gz] [R2.fastq.gz (optional)]

This script aligns fastq file(s) to a given genome.

OPTIONS:
   -b GENOME_BUILD		Genome build (e.g. mm10)

   -i STAR_GENOME_INDEX Path to the STAR genome index directory
   	  					Only required if there isn't already one for the passed build at /project/MDL_ChIPseq/process_sequencing_data/genome/star_index
   				 
   -n NTHREADS			Number of parallel threads

   -o OUTFILE			Output file name
EOF
}

# parse arguments (more info at https://ndench.github.io/bash/parsing-bash-flags)
nthreads=1
while getopts “:b:i:n:o:f:” OPTION; do
	case $OPTION in
		b) build=$OPTARG ;;
		i) star_index_dir=$(readlink -f $OPTARG) ;; # full path
		n) nthreads=$OPTARG ;;
		o) outfile=$(readlink -f $OPTARG) ;;
		?) usage; exit 1 ;;
 	esac
done
shift $(( OPTIND - 1 ))
fastq=$@ # read trailing arguments

# check if index exists
[[ ! -e $star_index_dir ]] && star_index_dir=/project/MDL_ChIPseq/process_sequencing_data/genome/star_index/$build # if star_index_dir was not passed, set it according to passed build
if [[ ! -e $star_index_dir ]]; then # if it still does not exist, exit and prompt the user to create index first
	echo "STAR genome not found. Either pass the folder with the STAR genome index files or create it from a fasta file using star_index.sh"
	exit 1
fi

# throw error if not all mandatory arguments are passed
if [ -z ${build+x} ] || [ -z ${star_index_dir+x} ] || [ -z ${nthreads+x} ] || [ -z ${outfile+x} ] || [ -z ${fastq+x} ]; then
	usage
	exit 1
fi

# define variables
[[ $nthreads -gt 20 ]] && nthreads=20 # cap number of parallel threads to 20 as STAR throws an error with more than 20 threads
outfile_prefix=${outfile%.*} # /path/to/file
outdir=${outfile%/*} # /path/to
sample_name=${outfile_prefix##*/} # file
log_dir=${outdir}/log
mkdir -p $log_dir
tmp=/scratch/local2/${sample_name}_tmp_${USER}
[[ -e $tmp ]] && rm -rf $tmp # clear tmp folder if it exists

# align reads
STAR --runMode alignReads \
	 --genomeLoad NoSharedMemory \
	 --genomeDir $star_index_dir \
	 --runThreadN $nthreads \
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
