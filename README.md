# process_seq_data

###
### Pipeline for processing ChIP-seq and ATAC-seq data from fastq to bigwig
###

Workflow:

1. Go to the location in your project folder where the data is to be stored.

2. Download the data table template and fill in the specifics of your data.
   For in-house data: https://bit.ly/38VOMrc
   For GEO data: https://bit.ly/...

3. Run process_data.sh: ./process_data.sh data_dir=. data_table=[seqcore_data.xlsx | geo_data.xlsx] nthreads=1
