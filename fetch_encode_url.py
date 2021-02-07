#! /usr/bin/env python

# Script for fetching the URL to the fastq of an ENCODE experiment
# Adapted from https://www.biostars.org/p/321633/

import numpy as np, pandas as pd, requests, os, sys
from pprint import pprint

def get(resource, url='https://www.encodeproject.org/{}/?format=json', headers={'accept': 'application/json'}):
    return requests.get(url.format(resource), headers=headers).json()

def format(file):
    return 'https://www.encodeproject.org' + file['href']

def get_exp(exp_acc, biol_rep, read, sequencing_type):
    response = get(os.path.join('experiments/', exp_acc))
    for file in response['files']:
        if (file['file_type'] == 'fastq'):
            if (sequencing_type=='paired_end'):
                if (int(file['replicate']['biological_replicate_number'])==biol_rep) & (int(file[sequencing_type])==read):
                    return format(file)
            elif (sequencing_type=='single_end') & (int(file['replicate']['biological_replicate_number'])==biol_rep):
                return format(file)

def main():
    if not len(sys.argv) == 5:
        msg = """Usage: python get_encode_url.py <Accession Number> <Biological Replicate> <Read> <Sequencing Type>\n
Accession Number\tMust be the one for the entire EXPERIMENT, not for the library, and not for the single fastq file.
Biological Replicate\tInteger number of biological replicate (1, 2, ...)
Read\t\t\tInteger number of single- / paired-end read (1 or 2)
Sequencing Type\t\tsingle-end or paired-end\n
Example: python get_encode_url.py ENCSR255XTC 1 2 paired-end
"""
        print(msg)
        sys.exit(0)
    _, exp_acc, biol_rep, read, sequencing_type = sys.argv
    sequencing_type = sequencing_type.replace('-','_')
    
    enc_url=get_exp(exp_acc, int(biol_rep), int(read), sequencing_type)
    print(enc_url)
    return
     
if __name__ == '__main__':
    main()
