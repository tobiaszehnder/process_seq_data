# process_seq_data

Pipeline for processing ChIP-seq and ATAC-seq data from fastq to bigwig for in-house data or from GEO

---

Workflow:

1. Go to the location in your project folder where the data is to be stored.

2. Download the data table template and fill in the specifics of your data.
   https://bit.ly/38VOMrc

3. Save the data table in comma-separated csv format

4. Run the pipeline: prun mdl process_data -d DATA_DIR -t DATA_TABLE -n NTHREADS [-o OVERWRITE]
   For more details on how to run the pipeline, run 'prun mdl process_data' without arguments

---

Notes:

You don't need to specify the exact file name of your fastq files.
This script will recognize R1 and R2 fastq files from paired-end sequencing and process them together.
Files from multiple sequencing runs lie in different folders and have to be specified individually.

---

zehnder [at] molgen.mpg.de
https://github.com/tobiaszehnder/process_seq_data
