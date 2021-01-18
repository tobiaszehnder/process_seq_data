# process_seq_data

Pipeline for processing ChIP-seq and ATAC-seq data from fastq to bigwig


---

Workflow:

1. Go to the location in your project folder where the data is to be stored.

2. Download the data table template and fill in the specifics of your data.
   For in-house data: https://bit.ly/38VOMrc
   For GEO data: https://bit.ly/...

3. Save the data table in comma-separated csv format

4. Run process_data.sh: ./process_data.sh data_dir=. data_table=data_table.csv nthreads=1


---

Notes:

You do not need to specify the exact file name of your fastq files.
This script will recognize R1 and R2 fastq files from paired-end sequencing and process them together.
Files from multiple sequencing runs lie in different folders and have to be specified individually.

---


zehnder [at] molgen.mpg.de