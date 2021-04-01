# process_seq_data

Pipeline for processing ChIP-seq and ATAC-seq data from fastq to bigwig for in-house data or from GEO

Scripts:

	process_data
	star_align
	star_index

---
 
process_data

	1. Go to the location in your project folder where the data is to be stored.

   	2. Download the data table template and fill in the specifics of your data.
	   https://bit.ly/38VOMrc

   	3. Save the data table in comma-separated csv format.

	4. Enter the command: prun mdl process_data -d DATA_DIR -t DATA_TABLE [optional arguments: -p NTHREADS -s STAR_PARAMS -o OVERWRITE -n -u -b -a]
   	   For more details on how to run the script, run 'prun mdl process_data' without arguments.


star_align

	This script will be used automatically in process_data for aligning ChIP-seq reads.


star_index

	This script can be used for creating new STAR genome index files from fasta in case the desired files do not yet exist at /project/MDL_ChIPseq/data/genome/star_index/.
	Enter the command: prun mdl star_index genome_build genome.fa nthreads

---

Notes

	You don't need to specify the exact file name of your fastq files.
	This script will recognize R1 and R2 fastq files from paired-end sequencing and process them together.
	Files from multiple sequencing runs lie in different folders and have to be specified individually.

---

zehnder [at] molgen.mpg.de
https://github.com/tobiaszehnder/process_seq_data
