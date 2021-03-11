# function for flattening a list of lists
flatten = lambda l: [item for sublist in l for item in sublist]

# function for retrieving SRX for a given SRR
get_srx = lambda srr: "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=%s&result=read_run&fields=study_accession,sample_accession,experiment_accession,run_accession,tax_id,scientific_name,fastq_ftp,submitted_ftp,sra_ftp&format=tsv&download=true" %srr

def parse_data_table(data_table):
    # load modules
    import numpy as np, pandas as pd, urllib.request, collections, glob, re, os
    
    # load data table and split it into dataframes for mpi / geo /encode data
    df = pd.read_csv(data_table)
    df.rename(columns={'seqcore_folder / SRR / ENC' : 'id'}, inplace=True)
    df.insert(0, 'species', df.loc[:,'build'].apply(lambda x: ''.join([i for i in x if not i.isdigit()])))
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

    if mpi.shape[0] > 0:
        # define bam and fastq sample names
        mpi.insert(0, 'source', 'MPI')
        mpi.insert(0, 'sample', mpi.loc[:, ['feature','tissue','stage','build','condition','biological_replicate','library_number']].apply(lambda x: '_'.join(x.dropna().values), axis=1))
        mpi.insert(0, 'fastq_sample', mpi.loc[:, ['feature','tissue','stage','species','condition','biological_replicate','library_number']].apply(lambda x: '_'.join(x.dropna().values), axis=1))
        mpi.insert(0, 'data_dir', mpi.loc[:,'species'].apply(lambda x: '/project/MDL_ChIPseq/data/epigenome/%s' %x))
        mpi.insert(0, 'reads', mpi.loc[:,'sequencing_type'].apply(lambda x: ['R1','R2'] if x=='paired-end' else ['R1']))
        # mpi.insert(0, 'seqcore_folder', mpi.loc[:,'id'])

        # link fastqs to fastq directory. rows from multiple sequencing runs are identified by having the same sample name and the same library number but different flow cell numbers
        # and different seqcore folders, so they must be linked to the fastq directory such that the snakefile can later find the files of both runs to merge them
        # therefore, leave the flow cell in the link name, and the snakefile will later find multiple files with different flow cell numbers
        for idx,row in mpi.iterrows():
            files = flatten([glob.glob('%s/*%s*%s_%s*fastq.gz' %tuple(row[['id','library_number','flow_cell']].tolist() + [r])) for r in row['reads']])
            directory = '%s/fastq/MPI/%s/' %(row['data_dir'], row['sequencing_type'])
            if not os.path.exists(directory):
                os.makedirs(directory)
            links = ['%s/%s_%s_%s.fastq.gz' %(directory, row['fastq_sample'], row['flow_cell'], r) for r in row['reads']]
            for i in range(len(files)):
                os.symlink(files[i], links[i] + '.tmp')
                os.rename(links[i] + '.tmp', links[i]) # that way, existing links are overwritten (symlink doesn't overwrite)
        mpi = mpi.loc[:,['sample', 'fastq_sample', 'source', 'experiment','sequencing_type', 'species', 'flow_cell']]
        mpi = mpi.groupby(['sample', 'fastq_sample', 'source', 'experiment','sequencing_type','species'])['flow_cell'].apply(','.join).reset_index() # join techn. repl.
    else:
        mpi = pd.DataFrame()
        
    # --------------------------
    ## GEO: load SRX to identify potential multiple sequencing runs and take them together in the sample name (e.g. sample_SRR1_SRR2)
    # --------------------------

    if geo.shape[0] > 0:
        # define sample strings
        geo.insert(0, 'source', 'GEO')
        geo.insert(0, 'sample', geo.loc[:, ['feature','tissue','stage','build','condition','biological_replicate']].apply(lambda x: '_'.join(x.dropna().values), axis=1))
        geo.insert(0, 'srx', geo.id.apply(lambda srr: str(urllib.request.urlopen(get_srx(srr)).read()).split('\\n')[-2].split('\\t')[2]))

        # add SRR to sample name (multiple in case of multiple sequencing runs)
        # rows from multiple sequencing runs are identified by having the same sample name and the same SRX and are merged
        grp_idx = geo.groupby(['sample','srx']).apply(lambda x: x.index.values)
        grp_val = geo.groupby(['sample','srx']).apply(lambda x: '_'.join(x['id'].values))
        new_sample_names = flatten([[geo.loc[x, 'sample'] + '_' + grp_val[i] for x in idx] for i,idx in enumerate(grp_idx.values)])
        geo.loc[:, 'sample'] = new_sample_names
        geo.insert(0, 'fastq_sample', geo.apply(lambda x: re.sub(x['build'], x['species'], x['sample']), axis=1))
        geo = geo.drop('id', axis=1).drop_duplicates()
        geo = geo.loc[:,['sample','fastq_sample','source','experiment','sequencing_type','species']]
    else:
        geo = pd.DataFrame()

    # --------------------------
    ## ENCODE: add experiment accession number to sample name. file accession numbers will be retrieved in a specific script (fetch_encode_url.py)
    # --------------------------
    if enc.shape[0] > 0:
        enc.insert(0, 'source', 'ENCODE')
        enc.insert(0, 'sample', enc.loc[:, ['feature','tissue','stage','build','condition','biological_replicate','id']].apply(lambda x: '_'.join(x.dropna().values), axis=1))
        enc.insert(0, 'fastq_sample', enc.loc[:, ['feature','tissue','stage','species','condition','biological_replicate', 'id']].apply(lambda x: '_'.join(x.dropna().values), axis=1))
        enc = enc.loc[:,['sample','fastq_sample','source','experiment','sequencing_type','species']]
    else:
        enc = pd.DataFrame()
    
    # --------------------------
    ## concat and return data (and replace all occurrences of dots)
    df = pd.concat([mpi, geo, enc], axis=0, sort=False).reset_index(drop=True).replace({'\.': ''}, regex=True)
    
    return df
