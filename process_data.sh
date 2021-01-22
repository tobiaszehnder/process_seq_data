#! /bin/bash

### Script to process ChIP-seq data from fastq to bigwig

# add program PATHs to environmental variable
export PATH=/home/zehnder/.local/bin/:$PATH

# # Define usage message
# usage="
# Usage: ./process_data.sh [--data_dir value] [--data_table value] [--nthreads value]\n
# -d\t--data_dir\tPath to your project data folder\n
# -t\t--data_table\tPath to your data table in comma-separated csv format\n
# \t\t\tDownload a template here:\n
# \t\t\tFor in-house data: https://bit.ly/38VOMrc\n
# \t\t\tFor GEO data: https://bit.ly/3sFS6OG\n
# -n\t--nthreads\tNumber of parallel threads\n
# -o\t--overwrite_file`s\tFlag for overwriting existing files.\n
# "

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

nthreads=1
overwrite=False
# Information about parsing bash arguments: https://ndench.github.io/bash/parsing-bash-flags
while getopts “:d:t:n:o” OPTION; do
	case $OPTION in
		d)
			data_dir=$OPTARG
			;;
		t)
			data_table=$OPTARG
			;;
		n)
			nthreads=$OPTARG
			;;
		o)
			overwrite=True
			;;
		?)
		usage
		exit 1
		;;
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
declare -A fastqs bams bigwigs experiments
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
		experiment=${arr[8]}
		experiments[$sample]=$experiment
		if [[ $sequencing_type == 'paired-end' ]]; then reads=( "R1" "R2" ); elif [[ $sequencing_type == 'single-end' ]]; then reads=( "R1" ); else echo "sequencing_type must be either 'single-end' or 'paired-end'" && exit 1; fi
		
		# create symbolic links for original fastq files
		files=()
		for read in "${reads[@]}"; do
			file=$seqcore_link/mpimg_$library_number*$read*fastq.gz
			link=${data_dir}/fastq/${feature}_${tissue}_${stage}_${build}_${flow_cell}_${read}.fastq.gz
			if [[ $overwrite==True ]] || [[ ! -e $link ]]; then
				ln -sf $file $link
			else
				echo "$link exists."
			fi
			# store the current link in an associative array with the name of the merged fastq as the key (in case of mutliple sequencing runs). Adding a comma allows to split them later.
			fastq=${data_dir}/fastq/${feature}_${tissue}_${stage}_${build}_${read}.fastq.gz
			fastqs[$fastq]+=$link,
		done
	done
} < $data_table

# merge fastq's of potential multiple sequencing runs
for fastq in "${!fastqs[@]}"; do
	if [[ $overwrite==True ]] || [[ ! -e $fastq ]]; then
		n=$(awk -F, '{print NF-1}' <<< "${fastqs[$fastq]}") # count the number of commas (i.e. the number of sequencing runs)
		if [[ $n == 1 ]]; then # one run: rename, remove flow-cell number
			[[ ! -f $fastq ]] && mv ${fastqs[$fastq]//$','/} $fastq
		elif [[ $n == 2 ]]; then # multiple runs: merge and remove flow-cell numbers
			echo -ne "\rMerging fastq files from multiple sequencing runs"
			[[ ! -f $fastq ]] && cat ${fastqs[$merged]//$','/$' '} > $fastq
		fi
	else
		echo "$fastq exists."
	# store the current fastq in an associative array with the bamfile as the key. Adding a comma allows to split potential paired-end fastqs later (R1 and R2).
	bam=$(sed -e 's/_R.\.fastq.gz/\.bam/g' -e 's/fastq/bam/g' <<< "$fastq")
	bams[$bam]+=$fastq,
done
echo ""

# Align data and produce bigwig tracks
for bam in "${!bams[@]}"; do
	echo ${bam//'\.bam'/':'}
	# align fastq to bam
	echo "Align reads"
	sample=$(sed -e 's/\.bam//g' <<< "$bam" | rev | cut -f1 -d/ | rev)
	build=$(echo $sample | cut -f4 -d_)
	if [[ $overwrite==True ]] || [[ ! -e $bam ]]; then
		star_align.sh $build $nthreads $bam ${bams[$bam]//','/' '} # the 'bams' array holds the associated fastqs (R1 / R2 for paired-end) at the entry with the bam name as key
	else
		echo "$bam exists."
	fi

	# sort bam and remove duplicates
	bam_sort=${bam//'\.bam'/'\.sort.bam'}
	bam_rmdup=${bam//'\.bam'/'\.sort.rmdup.bam'}
	echo "Sort alignment"
	if [[ $overwrite==True ]] || [[ ! -e $bam_sort ]]; then  
		samtools sort --threads $nthreads $bam > $bam_sort
	fi
	echo "Remove duplicates"
	if [[ $overwrite==True ]] || [[ ! -e $bam_rmdup ]]; then  
		samtools rmdup $bam_sort $bam_rmdup
	fi
	
	# index bam-file
	echo "Index bam-file"
	if [[ $overwrite==True ]] || [[ ! -e $bam.csi ]]; then
		samtools index -c $bam_rmdup # index bam-file: use -c for .csi format instead of .bai, allowing for chromosome sizes larger than 500 Mb (e.g. carPer2.1)
	fi

	# produce bigwig track
	echo "Produce bigwig track"
	bigwig=${echo $(sed -e 's/\.bam/\.cpm.bigwig' <<< "$bam")}
	if [[ exerpiments[$sample] == "ChIP-seq" ]]; then center_reads_flag="--center_reads"; else center_reads_flag=""; fi # center reads for ChIP-seq experiments, not for chromatin-accessiblity
	if [[ $overwrite==True ]] || [[ ! -e $bigwig ]]; then
		bamCoverage -binSize 10 --normalizeUsing RPKM $center_reads_flag --minMappingQuality 30 -p $nthreads -b $bam_rmdup -o $bigwig
	fi
done

### CLEANUP
# remove merged fastq's, only keep links to seqcore
# find $data_dir/fastq/ -type f -delete
# remove any intermediate bam files, only keep final sorted and rmduped bam
# find bam/ -type f -name *.bam ! -name '*.sort.rmdup.bam' -delete

echo Done
