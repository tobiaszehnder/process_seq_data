# function for flattening a list of lists
flatten = lambda l: [item for sublist in l for item in sublist]

# function for retrieving SRX for a given SRR
get_srx = lambda srr: "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=%s&result=read_run&fields=study_accession,sample_accession,experiment_accession,run_accession,tax_id,scientific_name,fastq_ftp,submitted_ftp,sra_ftp&format=tsv&download=true" %srr

def remove_mate_flags_function(path_infile, path_outfile):
    # function to remove mate flags from paired-end ATAC-seq bam files.
    # this is used as a preparation before creating the bigwig track so that bamCoverage treats the mates as single-end (because we're interested in the 5' cutting sites of the mates, not the read center)
    import pysam
    infile = pysam.AlignmentFile(str(path_infile), 'rb')
    outfile = pysam.AlignmentFile(str(path_outfile), 'wb', template=infile)
    for read in infile.fetch():
        read.flag = (read.flag & 3860) # 3860 are all flags which are not related to paired-end reads, i.e. a read flag will be either that or zero.
        outfile.write(read)
    infile.close()
    outfile.close()
    return

def check_and_standardize_data_table_values(df):
  import sys, os
  # feature
  df.loc[df.feature.str.lower().str.contains('atac'), 'feature'] = 'ATAC-seq'
  df.loc[df.feature.str.lower().str.contains('k27ac'), 'feature'] = 'H3K27ac'
  df.loc[df.feature.str.lower().str.contains('k4me1'), 'feature'] = 'H3K4me1'
  df.loc[df.feature.str.lower().str.contains('k4me3'), 'feature'] = 'H3K4me3'
  df.loc[df.feature.str.lower().str.contains('k27me3'), 'feature'] = 'H3K27me3'
  df.loc[df.feature.str.lower().str.contains('k9me3'), 'feature'] = 'H3K9me3'
  df.loc[df.feature.str.lower().str.contains('k36me3'), 'feature'] = 'H3K36me3'
  df.loc[df.feature.str.lower().str.contains('ctcf'), 'feature'] = 'CTCF'
  df.loc[df.feature.str.lower().str.contains('rad.*21'), 'feature'] = 'RAD21'
  df.loc[df.feature.str.lower().str.contains('input'), 'feature'] = 'Input'
  df['feature'] = df.feature.apply(lambda x: x[0].upper() + x[1:]) # always capitalize first letter
  df['biological_replicate'] = df.biological_replicate.str.replace('^[Rr][Ee][Pp][^0-9]*', 'Rep', regex=True) # Rep1, Rep2
  
  # seqcore folder / accession number
  if not df.id.apply(lambda x: x.startswith(('SRR','ENC')) or os.path.isdir(x)).all():
    raise ValueError('All entries in column "seqcore_folder / SRR / ENC" must be either a valid folder or an accession number starting with SRR or ENC')
    sys.exit(0)

  # sequencing type
  df.loc[df.sequencing_type.str.lower().str.contains('pair|pe'), 'sequencing_type'] = 'paired-end'
  df.loc[df.sequencing_type.str.lower().str.contains('single|se'), 'sequencing_type'] = 'single-end'
  if not df.sequencing_type.isin(['single-end','paired-end']).all():
    raise ValueError('All entries in column "sequencing_type" must be either "single-end" or "paired-end"')

  # experiment
  df.loc[df.experiment.str.lower().str.contains('chip.*seq', regex=True), 'experiment'] = 'ChIP-seq'
  df.loc[df.experiment.str.lower().str.contains('acc|atac', regex=True), 'experiment'] = 'chromatin-accessibility'
  if not df.experiment.isin(['ChIP-seq','chromatin-accessibility']).all():
    raise ValueError('All entries in column "experiment" must be either "ChIP-seq" or "chromatin-accessibility"')
    sys.exit(0)
    
  return df
  
def parse_data_table(data_table):
    # load modules
    import numpy as np, pandas as pd, urllib.request, collections, glob, re, sys, os
    
    # load data table and split it into dataframes for mpi / geo /encode data
    df = pd.read_csv(data_table, comment='#').replace({'\.': ''}, regex=True) # replace all occurrences of dots
    df.columns = ['id','feature','tissue','stage','build','condition','biological_replicate','sequencing_type','experiment','library_number','flow_cell']
    if df.loc[:,['id','feature','tissue','stage','build','condition','biological_replicate','sequencing_type','experiment']].isna().any().any():
      print('Missing values in data table. Only the columns `library_number` and `flow_cell` are allowed to be empty for non-MPI data entries.')
      sys.exit(0)
    df.insert(0, 'species', df.loc[:,'build'].apply(lambda x: ''.join([i for i in x if not i.isdigit()])))
    df = check_and_standardize_data_table_values(df)
    mpi_idx = df.loc[:,'id'].apply(lambda x: os.path.isdir(x))
    geo_idx = df.loc[:,'id'].apply(lambda x: x.startswith('SRR'))
    enc_idx = df.loc[:,'id'].apply(lambda x: x.startswith('ENC'))
    mpi = df.loc[mpi_idx]
    geo = df.loc[geo_idx]
    enc = df.loc[enc_idx]
    not_valid = df.loc[~(mpi_idx | geo_idx | enc_idx)]
    if not_valid.shape[0] > 0:
        print('The following rows do not have a valid entry at the first column (Path to seqcore folder, SRR or ENC accession number)')
        print(not_valid)

    # --------------------------
    ## MPI: create symlinks to fastq files
    # --------------------------

    ## MPI: write table with original seqcore file names and the associated sample names which the Snakefile can read
    if mpi.shape[0] > 0:
      # define bam and fastq sample names
      mpi.insert(0, 'source', 'MPI')
      mpi.insert(0, 'sample', mpi.loc[:, ['feature','tissue','stage','build','condition','biological_replicate','library_number']].apply(lambda x: '_'.join(x.values), axis=1))
      mpi.insert(0, 'fastq_sample', mpi.loc[:, ['feature','tissue','stage','species','condition','biological_replicate','library_number']].apply(lambda x: '_'.join(x.values), axis=1))
      mpi.insert(0, 'data_dir', mpi.loc[:,'species'].apply(lambda x: '/project/MDL_ChIPseq/data/epigenome/%s' %x))
      mpi.insert(0, 'mates', mpi.loc[:,'sequencing_type'].apply(lambda x: ['R1','R2'] if x=='paired-end' else ['R1']))

      # write links to original seqcore files to file
      file_dict = {}
      for idx,row in mpi.iterrows():
        for mate in row['mates']:
          file_dict['%s_%s_%s' %(row['fastq_sample'], row['flow_cell'], mate)] = glob.glob('%s/*%s*%s_%s*fastq.gz' %tuple(row[['id','library_number','flow_cell']].tolist() + [mate]))
      mpi_files_df = pd.DataFrame(file_dict, index=['original_file']).T
      mpi_files_df.to_csv(re.sub('.csv','.mpi_files',data_table))

      mpi = mpi.loc[:,['sample', 'fastq_sample', 'source', 'experiment','sequencing_type', 'species', 'build', 'flow_cell']]

      # join techn. repl.
      mpi = mpi.groupby(['sample', 'fastq_sample', 'source', 'experiment','sequencing_type','species','build'])['flow_cell'].apply(','.join).reset_index()
    else:
      mpi = pd.DataFrame()
    
    # --------------------------
    ## GEO: load SRX to identify potential multiple sequencing runs and take them together in the sample name (e.g. sample_SRR1_SRR2)
    # --------------------------

    if geo.shape[0] > 0:
        # define sample strings
        geo.insert(0, 'source', 'GEO')
        geo.insert(0, 'sample', geo.loc[:, ['feature','tissue','stage','build','condition','biological_replicate']].apply(lambda x: '_'.join(x.values), axis=1))
        geo.insert(0, 'srx', geo.id.apply(lambda srr: str(urllib.request.urlopen(get_srx(srr)).read()).split('\\n')[-2].split('\\t')[2]))

        # add SRR to sample name (multiple in case of multiple sequencing runs)
        # rows from multiple sequencing runs are identified by having the same sample name and the same SRX and are merged
        new_sample_names = pd.DataFrame(geo.groupby(['sample','srx']).apply(lambda x: '_'.join(x['id'])).reset_index().set_index('sample',drop=False).loc[:,['sample',0]].agg('_'.join, axis=1))
        geo_grouped = geo.groupby(['sample','srx']).apply(lambda x: x).set_index('sample').sort_values('sample')
        geo = pd.merge(new_sample_names, geo_grouped, on='sample').reset_index(drop=True)
        geo.rename(columns={0 : 'sample'}, inplace=True)
        geo.insert(1, 'fastq_sample', geo.apply(lambda x: re.sub(x['build'], x['species'], x['sample']), axis=1))
        geo = geo.loc[:,['sample','fastq_sample','source','experiment','sequencing_type','species','build']]
        geo = geo.drop_duplicates()
    else:
        geo = pd.DataFrame()

    # --------------------------
    ## ENCODE: add experiment accession number to sample name. file accession numbers will be retrieved in a specific script (fetch_encode_url.py)
    # --------------------------
    if enc.shape[0] > 0:
        enc.insert(0, 'source', 'ENCODE')
        enc.insert(0, 'sample', enc.loc[:, ['feature','tissue','stage','build','condition','biological_replicate','id']].apply(lambda x: '_'.join(x.values), axis=1))
        enc.insert(0, 'fastq_sample', enc.loc[:, ['feature','tissue','stage','species','condition','biological_replicate', 'id']].apply(lambda x: '_'.join(x.values), axis=1))
        enc = enc.loc[:,['sample','fastq_sample','source','experiment','sequencing_type','species','build']]
    else:
        enc = pd.DataFrame()
    
    # --------------------------
    ## concat data and replace all occurrences of dots
    df = pd.concat([mpi, geo, enc], axis=0, sort=False).reset_index(drop=True)
    
    return df
