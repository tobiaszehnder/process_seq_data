#! /bin/bash

### Script to process ChIP-seq data from fastq to bigwig

# . /project/wig-data/usr/profile

# parse arguments
for arg in "$@"
do
    key=$(echo $arg | cut -f1 -d=)
	val=$(echo $arg | cut -f2 -d=)   
    case "$key" in
		data_dir)   data_dir=${val} ;;
		data_table) data_table=${val} ;;
		nthreads)   nthreads=${val} ;;
		*)   
    esac    
done

usage_msg="Usage: ./process_data.sh data_dir=/path/to/your/project/data/folder data_table=/path/to/your/data/table nthreads=number_of_parallel_threads"
[ -z ${data_dir+x} ] || [ -z ${data_table+x} ] || [ -z ${nthreads+x} ] && echo $usage_msg && exit 1

# create folder structures
mkdir -p $data_dir/fastq $data_dir/bam $data_dir/bigwig
# shopt -s nullglob # prevents an error in case a glob does not match any name

# parse data table and create fastq links
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
			library_number=${arr[8]}
			flow_cell=${arr[9]}
		elif [[ ${header[0]} == "SRR_number" ]]; then # parse GEO data table
			SRR_number=${arr[0]}
		fi
		feature=${arr[1]}
		tissue=${arr[2]}
		stage=${arr[3]}
		build=${arr[4]}
		condition=${arr[5]}
		biological_replicate=${arr[6]}
		sequencing_type=${arr[7]}

		# create symbolic links for original fastq files
  	  	files=($seqcore_link/mpimg_$library# _number*R1*fastq.gz)
		[[ $sequencing_type == 'paired-end' ]] && files[1]=$seqcore_link/mpimg_$library_number*R2*fastq.gz # add R2 file for paired-end
		for j in ${!files[@]}; do # loop through array indices (0-based)
			link=${data_dir}/fastq/${feature}_${tissue}_${stage}_${build}_${flow_cell}_R$((j+1)).fastq.gz
			[[ ! -e $link ]] && ln -s $files[$i]} $link
		done
	done
} < $data_table

# merge fastq's of potential multiple sequencing runs

  
# process data
