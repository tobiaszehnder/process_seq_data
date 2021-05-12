# import modules
from pytools.persistent_dict import PersistentDict
import pandas as pd, os
from python_funcs import *

# parse arguments
project_dir = os.path.abspath(config['project_dir'].rstrip('/')) + '/' + os.path.basename(os.path.splitext(config['data_table'])[0])
unmapped_flag = "-u" if config['write_unaligned'] == True else ''
star_params = '' if not 'star_params' in config.keys() else config['star_params']
ucsc_folder = os.path.abspath(config['ucsc_folder'])

# define global variables
data_dir = '/project/MDL_ChIPseq/data/epigenome' if config['local'] == False else project_dir + '/data_local' #  + os.path.basename(os.path.splitext(config['data_table'])[0])
print('Reading data table')
data_df = parse_data_table(config['data_table']).set_index('sample', drop=False)
# data_df.to_csv('test.df')
data_df_fastq = data_df.set_index('fastq_sample', drop=False)
mates = {'single-end': ['R1'], 'paired-end': ['R1', 'R2']}
star_index_dir = '/project/MDL_ChIPseq/data/genome/star_index'
bowtie2_index_dir = '/project/MDL_ChIPseq/data/genome/bowtie2_index'

# read mpi_files table (if data_table involves MPI data)
mpi_files = re.sub('.csv','.mpi_files',config['data_table'])
if os.path.isfile(mpi_files):
    mpi_files_df = pd.read_csv(mpi_files, index_col=0)
    os.remove(mpi_files) # delete file after reading


# rules
# ---------------------------

def get_log_files(project_dir, df):
    logfiles = []
    for idx, row in df.iterrows():
        prefix = '%s/log/%s' %(project_dir, row['sample'])
        logfiles += [prefix + suffix for suffix in ('.full.bam.flagstat', '.rmdup.log', '.bamCoverage.log')]
        if row['experiment'] == 'chromatin-accessibility':
            logfiles += [prefix + '.full.bam.bowtie2.log']
        else:
            logfiles += [prefix + '.full.Log.final.out']
    return logfiles

def get_trackdb_targets(project_dir, ucsc_folder):
    if os.path.isdir(ucsc_folder):
        return ['%s/trackDb_%s.txt' %(project_dir, build) for build in set(data_df.build)]
    else:
        return []

def get_targets(wc):
    bw = expand('%s/bigwig/{sample}.cpm.bw' %project_dir, sample=data_df['sample'])
    bam = expand('%s/bam/{sample}.rmdup.bam' %project_dir, sample=data_df['sample'])
    csi = expand('%s/bam/{sample}.rmdup.bam.csi' %project_dir, sample=data_df['sample'])
    logs = get_log_files(project_dir, data_df.loc[:,['sample','experiment']])
    trackdb = get_trackdb_targets(project_dir, ucsc_folder)
    targets = bw + bam + csi + logs + trackdb
    return targets
    
rule all:
    input:
        get_targets

rule link_seqcore:
    # here, sample includes flow cell and mate.
    input:
        lambda wc: mpi_files_df.loc[wc['sample'], 'original_file']
    output:
        '{data_dir}/fastq/MPI/{sequencing_type}/{sample}.fastq.gz'
    shell:
        'ln -sf {input} {output}'
        
def get_geo_fastqs(wc):
    # This function takes a fastq sample name with all trailing SRRs and returns a list of individual fastqs
    sample = '_'.join([wc[x] for x in ['feature', 'tissue', 'stage', 'build', 'condition', 'biological_replicate', 'id']])
    sample_without_id = '_'.join([wc[x] for x in ['feature', 'tissue', 'stage', 'build', 'condition', 'biological_replicate']])
    srrs = wc['id'].split('_')
    if data_df.loc[sample, 'sequencing_type'] == 'single-end':
        fastqs = ['%s/fastq/GEO/single-end/%s_%s_R1.fastq.gz' %(wc['data_dir'], sample_without_id, srr) for srr in srrs]
    elif data_df.loc[sample, 'sequencing_type'] == 'paired-end':
        fastqs = flatten([['%s/fastq/GEO/paired-end/%s_%s_%s.fastq.gz' %(wc['data_dir'], sample_without_id, srr, mate) for mate in ['R1','R2']] for srr in srrs])
    return fastqs

rule download_geo:
    output:
        temp('{prefix}/{sample}_{srr}.sra')
    shell:
        'prefetch -p -f yes -o {output} {wildcards.srr}'

rule fasterq_single_end:
    input:
        '{data_dir}/fastq/GEO/single-end/{sample}_SRR{srr_number}.sra'
    output:
        temp('{data_dir}/fastq/GEO/single-end/{sample}_SRR{srr_number}_1.fastq')
    threads: min(workflow.cores, 10)
    params:
        ncbi_tmp = '/scratch/local2/$USER/ncbi-tmp',
        output = '{data_dir}/fastq/GEO/single-end/{sample}_SRR{srr_number}.fastq',
        wrong_paired_end_mate2 = '{data_dir}/fastq/GEO/single-end/{sample}_SRR{srr_number}_2.fastq'
    shell:
        '''
        mkdir -p {params.ncbi_tmp}
        fasterq-dump -t {params.ncbi_tmp} --split-files -f -e {threads} -p -o {params.output} {input}
        if [[ -e {params.wrong_paired_end_mate2} ]]; then echo "Skipping SRR{wildcards.srr_number} (paired-end, but stated as single-end)"; rm {input} {output} {params.wrong_paired_end_mate2}; exit 0; fi
        mv {params.output} {output}
        '''
rule fasterq_paired_end:
    input:
        '{data_dir}/fastq/GEO/paired-end/{sample}_SRR{srr_number}.sra'
    output:
        R1 = temp('{data_dir}/fastq/GEO/paired-end/{sample}_SRR{srr_number}_1.fastq'),
        R2 = temp('{data_dir}/fastq/GEO/paired-end/{sample}_SRR{srr_number}_2.fastq')
    threads: min(workflow.cores, 10)
    params:
        ncbi_tmp = '/scratch/local2/$USER/ncbi-tmp',
        output = '{data_dir}/fastq/GEO/paired-end/{sample}_SRR{srr_number}.fastq'
    shell:
        """
        mkdir -p {params.ncbi_tmp}
        fasterq-dump -t {params.ncbi_tmp} --split-files -f -e {threads} -p -o {params.output} {input}
        """
        
rule zip_geo_fastq:
    input:
        '{seqtype_dir}/{sample}_SRR{srr_number}_{mate_number}.fastq'
    output:
        temp('{seqtype_dir}/{sample}_SRR{srr_number,[0-9]+}_{mate_number,[0-9]+}.fastq.gz')
    shell:
        'gzip -1 {input}'

rule rename_geo_fastq:
    input:
        '{seqtype_dir}/{sample}_SRR{srr_number}_{mate_number}.fastq.gz'
    output:
        temp('{seqtype_dir}/{sample}_SRR{srr_number,[0-9]+}_R{mate_number,[0-9]+}.fastq.gz')
    shell:
        'mv {input} {output}'

rule merge_geo_fastqs:
    input:
        lambda wc: ['{data_dir}/fastq/GEO/{sequencing_type}/{feature}_{tissue}_{stage}_{build}_{condition}_{biological_replicate}_%s_{mate}.fastq.gz' %id for id in wc['ids'].split('_')]
    output:
        temp('{data_dir}/fastq/GEO/{sequencing_type}/{feature,((?!_).)*}_{tissue,((?!_).)*}_{stage,((?!_).)*}_{build,((?!_).)*}_{condition,((?!_).)*}_{biological_replicate,((?!_).)*}_{ids}_merged_{mate}.fastq.gz')
    shell:
        'cat {input} > {output}'

rule merge_mpi_fastqs:
    # flow cell info is stored in data_df
    input:
        lambda wc: [ancient('%s/fastq/MPI/%s/%s_%s_%s.fastq.gz') %(wc['data_dir'], wc['sequencing_type'], wc['sample'], fc, wc['mate']) for fc in data_df_fastq.loc[wc['sample'], 'flow_cell'].split(',')]
    output:
        temp('{data_dir}/fastq/MPI/{sequencing_type}/{sample}_merged_{mate,\w\d}.fastq.gz')
    shell:
        'cat {input} > {output}'
        
rule download_encode_single_end:
    output:
        temp('{data_dir}/fastq/ENCODE/single-end/{feature}_{tissue}_{stage}_{build}_{condition}_{biological_replicate}_{id}_R1.fastq.gz')
    shell:
        """
        url=$(prun mdl fetch_encode_url.py {wildcards.id} $(sed -e 's/[^0-9]//g' <<< "{wildcards.biological_replicate}") 1 single-end)
        curl -o {output} -O -L $url
        """

rule download_encode_paired_end:
    output:
        R1 = temp('{data_dir}/fastq/ENCODE/paired-end/{feature}_{tissue}_{stage}_{build}_{condition}_{biological_replicate}_{id}_R1.fastq.gz'),
        R2 = temp('{data_dir}/fastq/ENCODE/paired-end/{feature}_{tissue}_{stage}_{build}_{condition}_{biological_replicate}_{id}_R2.fastq.gz')
    shell:
        """
        url_1=$(prun mdl fetch_encode_url.py {wildcards.id} $(sed -e 's/[^0-9]//g' <<< "{wildcards.biological_replicate}") 1 paired-end)
        url_2=$(prun mdl fetch_encode_url.py {wildcards.id} $(sed -e 's/[^0-9]//g' <<< "{wildcards.biological_replicate}") 2 paired-end)
        curl -o {output.R1} -O -L $url_1
        curl -o {output.R2} -O -L $url_2
        """

rule trim_adapters:
    input:
        # ancient: ignore time stamp of fastqs because in case of MPI data, potentially existing links are updated every time a data table is parsed, thereby evoking the execution of all downstream rules.
        R1 = ancient('{data_dir}/fastq/{source}/paired-end/{sample}_R1.fastq.gz'),
        R2 = ancient('{data_dir}/fastq/{source}/paired-end/{sample}_R2.fastq.gz')
    output:
        R1 = temp('{data_dir}/fastq/{source}/paired-end/{sample,ATAC((?!merged).)*}_R1_trim.fastq.gz'),
        R2 = temp('{data_dir}/fastq/{source}/paired-end/{sample,ATAC((?!merged).)*}_R2_trim.fastq.gz')
    log:
        '{data_dir}/bam/{source}/paired-end/{sample}_fastq_cutadapt.log'
    params:
        adapt1 = 'CTGTCTCTTATACACATCTCCGAGCCCACGAGAC',
        adapt2 = 'CTGTCTCTTATACACATCTGACGCTGCCGACGA'
    threads: min(workflow.cores, 10)
    shell:
        'cutadapt -j {threads} -a {params.adapt1} -A {params.adapt2} -o {output.R1} -p {output.R2} {input.R1} {input.R2} &> {log}'

rule star_index:
    input:
        ancient('/project/MDL_ChIPseq/data/genome/fasta/{build}.fa') # without ancient, this rule will be executed if fasta is newer than star_index_dir
    output:
        directory('%s/{build}' %star_index_dir)
    threads: min(workflow.cores, 10)
    shell:
        'prun mdl star_index {wildcards.build} {input} {threads}'

rule bowtie2_index:
    input:
        fasta = ancient('/project/MDL_ChIPseq/data/genome/fasta/{build}.fa'),
        sizes = ancient('/project/MDL_ChIPseq/data/genome/assembly/{build}.sizes')
    output:
        directory('%s/{build}' %bowtie2_index_dir)
    log:
        '%s/{build}/{build}.bowtie2-build.log' %bowtie2_index_dir
    threads: min(workflow.cores, 10)
    shell:
        '''
        # if [[ $(cat {input.sizes} | cut -f 2 | paste -sd+ | bc) < 4000000000 ]]; then large_index_flag="--large-index"; else large_index_flag=""; fi # only needed for forced large index for genomes < 4 Bbp
        mkdir -p {output}
        bowtie2-build --threads {threads} {input.fasta} {output}/{wildcards.build} &> {log}
        '''

def get_fastqs(wc):
    # This function creates input files for the rule `align` with source-specific folders (MPI/GEO/ENCODE).
    # In addition, it adds the suffix '_merged' for experiments with multiple replicates, which then evokes the respective rules to merge the individual fastqs.
    sample = '_'.join([wc[x] for x in ['feature', 'tissue', 'stage', 'build', 'condition', 'biological_replicate', 'id']])
    species = data_df.loc[sample, 'species']
    fastq_sample = re.sub(wc['build'], species, sample)
    if wc['source'] == 'MPI':
        suffix = '_merged' if len(data_df.loc[sample, 'flow_cell'].split(',')) > 1 else '_' + data_df.loc[sample, 'flow_cell']
    else:
         suffix = '_merged' if (len(wc['id'].split('_')) > 1) else ''
    if data_df.loc[sample, 'sequencing_type'] == 'single-end':
        fastqs = ['%s/fastq/%s/single-end/%s%s_R1.fastq.gz' %(wc['data_dir'], wc['source'], fastq_sample, suffix)]
    elif data_df.loc[sample, 'sequencing_type'] == 'paired-end':
        fastqs = ['%s/fastq/%s/paired-end/%s%s_%s.fastq.gz' %(wc['data_dir'], wc['source'], fastq_sample, suffix, mate) for mate in ['R1','R2']]
    if data_df.loc[sample, 'experiment'] == 'chromatin-accessibility' and data_df.loc[sample, 'sequencing_type'] == 'paired-end':
        fastqs = [re.sub('.fastq.gz', '_trim.fastq.gz', fastq) for fastq in fastqs]
    # ancient: ignore time stamp of fastqs because in case of MPI data, potentially existing links are updated every time a data table is parsed, thereby evoking the execution of all downstream rules.
    return [ancient(fastq) for fastq in fastqs]
        
rule align_star:
    input:
        index = '%s/{build}' %star_index_dir,
        fastqs = get_fastqs
    output:
        temp('{data_dir}/bam/{source}/{sequencing_type}/{feature,((?!ATAC|_).)*}_{tissue,((?!_).)*}_{stage,((?!_).)*}_{build,((?!_).)*}_{condition,((?!_).)*}_{biological_replicate,((?!_).)*}_{id}.full.bam')
    log:
        '{data_dir}/bam/{source}/{sequencing_type}/log/{feature}_{tissue}_{stage}_{build}_{condition}_{biological_replicate}_{id}.full.Log.final.out'
    params:
        star_params = "-s '%s'" %star_params if not star_params == '' else ''
    threads: min(workflow.cores, 10)
    shell:
        "prun mdl star_align {unmapped_flag} -b {wildcards.build} -i {input.index} -n {threads} {params.star_params} -o {output} {input.fastqs}"

rule align_bowtie2_single_end:
    # only for ATAC
    input:
        index = '%s/{build}' %bowtie2_index_dir,
        fastq = get_fastqs
    output:
        temp("{data_dir}/bam/{source}/single-end/{feature,ATAC.*}_{tissue}_{stage}_{build}_{condition}_{biological_replicate}_{id}.full.bam")
    log:
        '{data_dir}/bam/{source}/single-end/log/{feature,((?!_).)*}_{tissue,((?!_).)*}_{stage,((?!_).)*}_{build,((?!_).)*}_{condition,((?!_).)*}_{biological_replicate,((?!_).)*}_{id}.full.bam.bowtie2.log'
    threads: min(workflow.cores, 10)
    shell:
        'bowtie2 --mm -x {input.index}/{wildcards.build} --threads {threads} -U {input.fastq} 2> {log} | samtools view -Su /dev/stdin | samtools sort > {output}'

rule align_bowtie2_paired_end:
    # only for ATAC
    input:
        index = '%s/{build}' %bowtie2_index_dir,
        fastqs = get_fastqs
    output:
        temp("{data_dir}/bam/{source}/paired-end/{feature,ATAC((?!_).)*}_{tissue,((?!_).)*}_{stage,((?!_).)*}_{build,((?!_).)*}_{condition,((?!_).)*}_{biological_replicate,((?!_).)*}_{id}.full.bam")
    log:
        '{data_dir}/bam/{source}/paired-end/log/{feature}_{tissue}_{stage}_{build}_{condition}_{biological_replicate}_{id}.full.bam.bowtie2.log'
    threads: min(workflow.cores, 10)
    shell:
        'bowtie2 --mm -x {input.index}/{wildcards.build} --threads {threads} -1 {input.fastqs[0]} -2 {input.fastqs[1]} 2> {log} | samtools view -Su /dev/stdin | samtools sort > {output}'

rule flagstat:
    input:
        ancient('{data_dir}/bam/{source}/{sequencing_type}/{sample}.full.bam')
    output:
        '{data_dir}/bam/{source}/{sequencing_type,((?!log).)*}/log/{sample,((?!log).)*}.full.bam.flagstat'
    threads: min(workflow.cores, 10)
    shell:
        'samtools flagstat --threads {threads} {input} > {output}'
        
rule nsort:
    # only for paired-end because fixmate needs n-sorted bam
    input:
        bam = "{data_dir}/bam/{source}/paired-end/{sample}.full.bam",
        flagstat = "{data_dir}/bam/{source}/paired-end/log/{sample}.full.bam.flagstat"
    output:
        temp("{data_dir}/bam/{source}/paired-end/{sample}.nsort.bam")
    threads: min(workflow.cores, 10)
    shell:
        "samtools sort -n --threads {threads} -o {output} {input.bam}"

rule fixmate:
    input:
        "{data_dir}/bam/{source}/paired-end/{sample}.nsort.bam"
    output:
        temp("{data_dir}/bam/{source}/paired-end/{sample}.fixmate.bam")
    threads: min(workflow.cores, 10)
    shell:
        "samtools fixmate -m --threads {threads} {input} {output}"

rule csort:
    input:
        "{data_dir}/bam/{source}/paired-end/{sample}.fixmate.bam"
    output:
        temp("{data_dir}/bam/{source}/paired-end/{sample}.csort.bam")
    threads: min(workflow.cores, 10)
    shell:
        'samtools sort --threads {threads} -o {output} {input}'

rule rename_bam_single_end:
    # single-end full.bam files can be directly renamed to .csort.bam as they don't need to be fixmated and the output of STAR is already c-sorted.
    input:
        bam='{data_dir}/bam/{source}/single-end/{sample}.full.bam',
        flagstat="{data_dir}/bam/{source}/single-end/log/{sample}.full.bam.flagstat"
    output:
        temp('{data_dir}/bam/{source}/single-end/{sample}.csort.bam')
    shell:
        'mv {input.bam} {output}'
        
rule rmdup:
    input:
        "{data_dir}/bam/{source}/{sequencing_type}/{sample}.csort.bam"
    output:
        "{data_dir}/bam/{source}/{sequencing_type}/{sample}.rmdup.bam"
    log:
        "{data_dir}/bam/{source}/{sequencing_type}/log/{sample}.rmdup.log"
    threads: min(workflow.cores, 10)
    shell:
        "samtools markdup -s -r --threads {threads} {input} {output} 2> {log}"
        
rule bam_index:
    input:
        '{data_dir}/bam/{source}/{sequencing_type}/{sample}.rmdup.bam'
    output:
        '{data_dir}/bam/{source}/{sequencing_type}/{sample,((?!noMateFlags).)*}.rmdup.bam.csi'
    shell:
        "samtools index -c {input}"

rule bam_index_noMateFlag:
    # separate rule because this output is temporary
    input:
        '{data_dir}/bam/{source}/{sequencing_type}/{sample}.rmdup.noMateFlags.bam'
    output:
        temp('{data_dir}/bam/{source}/{sequencing_type}/{sample}.rmdup.noMateFlags.bam.csi')
    shell:
        "samtools index -c {input}"

rule remove_mate_flags:
    # remove mate flags for bigwig track creation for paired-end ATAC-seq (single-end mode)
    input:
        bam = '{source_dir}/paired-end/{sample}.rmdup.bam',
        idx = '{source_dir}/paired-end/{sample}.rmdup.bam.csi'
    output:
        temp('{source_dir}/paired-end/{sample,ATAC.*}.rmdup.noMateFlags.bam')
    run:
        remove_mate_flags_function(input.bam, output)

def get_bam(wc, index=False):
    suffix = '.csi' if index else ''
    if wc['sequencing_type'] == 'paired-end' and wc['sample'].startswith('ATAC'):
        return '%s/bam/%s/%s/%s.rmdup.noMateFlags.bam%s' %(wc['data_dir'], wc['source'], wc['sequencing_type'], wc['sample'], suffix)
    else:
        return '%s/bam/%s/%s/%s.rmdup.bam%s' %(wc['data_dir'], wc['source'], wc['sequencing_type'], wc['sample'], suffix)

rule bw:
    input:
        bam = get_bam,
        idx = lambda wc: get_bam(wc, index=True)
    output:
        "{data_dir}/bigwig/{source}/{sequencing_type}/{sample}.cpm.bw"
    log:
        "{data_dir}/bam/{source}/{sequencing_type}/log/{sample}.bamCoverage.log"
    params:
        center_reads_flag = lambda wc: '--centerReads' if data_df.loc[wc['sample'], 'experiment'] == 'ChIP-seq' else ''
    threads: min(workflow.cores, 10)
    shell:
        "bamCoverage --binSize 10 --normalizeUsing CPM {params.center_reads_flag} --minMappingQuality 30 -p {threads} -b {input.bam} -o {output} > {log}"

rule link_final_files:
    input:
        lambda wc: ancient('%s/%s/{filetype_dir}/%s/%s/{sample}.{ext1}.{ext2}' %(data_dir, data_df.loc[wc['sample'], 'species'], data_df.loc[wc['sample'], 'source'], data_df.loc[wc['sample'], 'sequencing_type']))
    output:
        '%s/{filetype_dir,((?!log).)*}/{sample,(?!log)[\w\-]+}.{ext1}.{ext2}' %project_dir
    shell:
        'ln -sf {input} {output}'

rule link_log_files:
    input:
        lambda wc: ancient('%s/%s/bam/%s/%s/log/{sample}.{ext1}.{ext2}' %(data_dir, data_df.loc[wc['sample'], 'species'], data_df.loc[wc['sample'], 'source'], data_df.loc[wc['sample'], 'sequencing_type']))
    output:
        '%s/log/{sample,[\w\-]+}.{ext1}.{ext2}' %project_dir
    shell:
        'ln -sf {input} {output}'

rule create_ucsc_trackdb:
    input:
        ucsc_folder = ucsc_folder,
        bigwigs = lambda wc: expand('%s/bigwig/{sample}.cpm.bw' %project_dir, sample=data_df.loc[data_df.build==wc.build, 'sample'])
    output:
        trackdb = '{project_dir}/trackDb_{build}.txt'
    shell:
        '''
        prun mdl create_ucsc_trackdb.py {output} {input.ucsc_folder} {input.bigwigs}
        cp {input.bigwigs} {input.ucsc_folder}
        '''
