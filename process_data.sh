#! /bin/bash

### Script to process ChIP-seq data from fastq to bigwig

# . /project/wig-data/usr/profile

# parse arguments
usage_msg="Usage: ./process_data.sh data_dir=/path/to/your/project/data/folder data_table=/path/to/your/data/table nthreads=10"
if [ -z ${data_dir+x} ] || [ -z ${data_table+x} ] || [ -z ${nthreads+x} ]; then
	echo usage_msg && exit 1
fi

# parse data table

# process data
