#! /bin/bash

### Script to process ChIP-seq data from fastq to bigwig

set -e # exit when any command fails
export PATH=/home/zehnder/.local/bin/:$PATH # add program PATHs to environmental variable

# ---------------------------------------------------
### Functions
# ---------------------------------------------------

usage()
{
cat << EOF
usage: $0 -d DATA_DIR -t DATA_TABLE -n NTHREADS [-o]

This script processes sequencing data from ChIP-seq and chromatin accessibility assays from fastq to bigwig.

OPTIONS:
   -d DATA_DIR   Path to your project data folder

   -t DATA_TABLE Path to your data table in comma-separated csv format
   	  	 Download a template here:
		 For in-house data: https://bit.ly/38VOMrc
		 For GEO data: https://bit.ly/3sFS6OG

   -n NTHREADS   Number of parallel threads

   -o OVERWRITE	 Flag for overwriting existing files
EOF
}

file_exists() {
	file=$(echo $1 | rev | cut -d'/' -f1,2 | rev) # trimmed path with only the last folder to avoid overcrowded print
    echo "$file  exists and will not be overwritten. (Set the -o option for overwriting existing files)"
}

# ---------------------------------------------------

# parse arguments (more info at https://ndench.github.io/bash/parsing-bash-flags)
nthreads=1
overwrite=False
while getopts “:d:t:n:o” OPTION; do
	case $OPTION in
		d) data_dir=$(readlink -f $OPTARG) ;; # full path
		t) data_table=$OPTARG ;;
		n) nthreads=$OPTARG ;;
		o) overwrite=True ;;
		?) usage; exit 1 ;;
 	esac
done

# throw error if not all mandatory arguments are passed
if [ -z ${data_dir+x} ] || [ -z ${data_table+x} ] || [ -z ${nthreads+x} ] || [ -z ${overwrite+x} ]; then
	usage
	exit 1
fi

# create folder structures
mkdir -p $data_dir/fastq $data_dir/bam $data_dir/bigwig

# parse data table and create fastq links
declare -A fastqs bams bigwigs experiments sequencing_types
{
	IFS="," read -r row || [ -n "$row" ] # read the header
	header=(${row//,/ }) # convert header string to array
	while IFS="," read -r row || [ -n "$row" ]; do # loop through rows
		# parse row from table and assign its values to variables
		row=${row//$'\r'/} # remove carriage returns
		row=${row//$'\n'/} # and newlines from row string
		arr=(${row//,/ }) # convert row string to array
		if [[ ${header[0]} == "seqcore_link" ]]; then # parse in-house data table
			seqcore_link=${arr[0]}
			library_number=${arr[9]}
			flow_cell=${arr[10]}
		elif [[ ${header[0]} == "SRR_number" ]]; then # parse GEO data table
			SRR_number=${arr[0]}
		fi
		feature=${arr[1]}
		tissue=${arr[2]}
		stage=${arr[3]}
		build=${arr[4]}
		sample=${feature}_${tissue}_${stage}_${build}
		condition=${arr[5]}
		biological_replicate=${arr[6]}
		sequencing_type=${arr[7]}
		sequencing_types[$sample]=$sequencing_type
		experiments[$sample]=${arr[8]}
		if [[ $sequencing_type == 'paired-end' ]]; then reads=( "R1" "R2" ); elif [[ $sequencing_type == 'single-end' ]]; then reads=( "R1" ); else echo "sequencing_type must be either 'single-end' or 'paired-end'" && exit 1; fi
		
		# create symbolic links for original fastq files
		files=()
		for read in "${reads[@]}"; do
			file=$seqcore_link/mpimg_${library_number}*${flow_cell}_${read}*fastq.gz
			link=${data_dir}/fastq/${feature}_${tissue}_${stage}_${build}_${flow_cell}_${read}.fastq.gz
			if [[ $overwrite == True ]] || [[ ! -e $link ]]; then
				ln -sf $file $link
			else
				file_exists $link
			fi
			# store the current link in an associative array with the name of the merged fastq as the key (in case of mutliple sequencing runs). Adding a comma allows to split them later.
			fastq=${data_dir}/fastq/${feature}_${tissue}_${stage}_${build}_${read}.fastq.gz
			fastqs[$fastq]+=$link,
		done
	done
} < $data_table

# merge fastq's of potential multiple sequencing runs
for fastq in "${!fastqs[@]}"; do
	if [[ $overwrite == True ]] || [[ ! -e $fastq ]]; then
		n=$(awk -F, '{print NF-1}' <<< "${fastqs[$fastq]}") # count the number of commas (i.e. the number of sequencing runs)
		if [[ $n == 1 ]]; then # one run: rename, remove flow-cell number
			mv ${fastqs[$fastq]//$','/} $fastq
		elif [[ $n == 2 ]]; then # multiple runs: merge and remove flow-cell numbers
			echo -ne "\rMerging fastq files from multiple sequencing runs"
			cat ${fastqs[$fastq]//$','/$' '} > $fastq
		fi
	else
		file_exists $fastq
	fi
	# store the current fastq in an associative array with the bamfile as the key. Adding a comma allows to split potential paired-end fastqs later (R1 and R2).
	bam=$(sed -e 's/_R.\.fastq.gz/\.bam/g' -e 's/fastq/bam/g' <<< "$fastq")
	bams[$bam]+=$fastq,
done
echo ""

# Align data and produce bigwig tracks
for bam in "${!bams[@]}"; do
	sample=$(sed -e 's/\.bam//g' <<< "$bam" | rev | cut -f1 -d/ | rev)
	build=$(echo $sample | cut -f4 -d_)
	echo -e "\n$sample:"

	# align fastq to bam
	echo "Align reads"
	if [[ $overwrite == True ]] || [[ ! -e $bam ]]; then
		# check if STAR genome index exists
		star_index_dir=/project/MDL_ChIPseq/process_sequencing_data/genome/star_index/$build # if star_index_dir was not passed, set it according to passed build
		if [[ ! -e $star_index_dir ]]; then # if it still does not exist, exit and prompt the user to create index first
			echo "STAR genome not found. Create it from a fasta file using star_index.sh"
			exit 1
		fi
		# align
		star_align.sh -b $build -n $nthreads -o $bam ${bams[$bam]//','/' '} # the 'bams' array holds the associated fastqs (R1 / R2 for paired-end) at the entry with the bam name as key
	else
		file_exists $bam
	fi

	# fill in mate coordinates for paired-end. this requires the bam to be sorted by name. later, markdup again requires the bam to be sorted by coordinate, which we need two steps of sorting here.
	if [[ ${sequencing_types[$sample]} == "paired-end" ]]; then
		# sort bam by name
		echo "Sort alignment by name"
		bam_nsort=${data_dir}/bam/${sample}.nsort.bam
		if [[ $overwrite == True ]] || [[ ! -e $bam_nsort ]]; then  
			samtools sort -n --threads $nthreads -o $bam_nsort $bam # -n flag for sorting by name
		else
			file_exists $bam_nsort
		fi

		# fixmate
		echo "Fill in mate coordinates"
		bam_fixmate=${data_dir}/bam/${sample}.fixmate.bam
		if [[ $overwrite == True ]] || [[ ! -e $bam_fixmate ]]; then  
			samtools fixmate -m --threads $nthreads $bam_nsort $bam_fixmate
		else
			file_exists $bam_fixmate
		fi

		# sort bam by coordinate
		echo "Sort alignment by coordinate"
		bam_csort=${data_dir}/bam/${sample}.csort.bam
		if [[ $overwrite == True ]] || [[ ! -e $bam_csort ]]; then  
			samtools sort --threads $nthreads -o $bam_csort $bam_fixmate
		else
			file_exists $bam_csort
		fi
		bam=$bam_csort
	fi
	
   # remove duplicates
	echo "Remove duplicates"
	bam_rmdup=${data_dir}/bam/${sample}.rmdup.bam
	if [[ $overwrite == True ]] || [[ ! -e $bam_rmdup ]]; then
		samtools markdup -r --threads $nthreads $bam $bam_rmdup
	else
		file_exists $bam_rmdup
	fi
	
	# index bam
	echo "Index bam-file"
	if [[ $overwrite == True ]] || [[ ! -e $bam_rmdup.csi ]]; then
		samtools index -c $bam_rmdup # index bam-file: use -c for .csi format instead of .bai, allowing for chromosome sizes larger than 500 Mb (e.g. carPer2.1)
	else
		file_exists ${bam}.csi
	fi

	# produce bigwig track
	echo "Produce bigwig track"
	bigwig=${data_dir}/bigwig/${sample}.rpkm.bw
 	if [[ exerpiments[$sample] == "ChIP-seq" ]]; then center_reads_flag="--center_reads"; else center_reads_flag=""; fi # center reads for ChIP-seq experiments, not for chromatin-accessiblity
	if [[ $overwrite == True ]] || [[ ! -e $bigwig ]]; then
		bamCoverage --binSize 10 --normalizeUsing RPKM $center_reads_flag --minMappingQuality 30 -p $nthreads -b $bam_rmdup -o $bigwig
	else
		file_exists $bigwig
	fi
done

# remove loaded genome from shared memory
# STAR=/scratch/ngsvin/bin/mapping/STAR/STAR_2.6.1d/bin/Linux_x86_64_static/STAR
STAR=STAR
$STAR --genomeDir $star_index_dir --genomeLoad Remove --outFileNamePrefix ${data_dir}/bam/log/removeGenomeLoad_${build}_

### CLEANUP
find $data_dir/fastq/ -type f -delete # remove merged fastq's, only keep links to seqcore
find bam/ -type f -name *.bam ! -name '*.rmdup.bam' -delete # remove any intermediate bam files, only keep final sorted and rmduped bam

echo Done
