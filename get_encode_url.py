#! /usr/bin/env python

# Script for fetching the URL to the fastq of an ENCODE experiment
# Adapted from https://www.biostars.org/p/321633/

import numpy as np, pandas as pd, requests, os, sys
from pprint import pprint

def get(resource, url='https://www.encodeproject.org/{}/?format=json', headers={'accept': 'application/json'}):
    return requests.get(url.format(resource), headers=headers).json()

def format(file):
    return 'https://www.encodeproject.org' + file['href']
  
def get_exp(exp_acc, biol_rep, read):
    response = get(os.path.join('experiments/', exp_acc))
    controls = set()
    for file in response['files']:
        if (file['file_type'] == 'fastq'):
            if (int(file['replicate']['biological_replicate_number'])==biol_rep) & (int(file['paired_end'])==read):
                return format(file)

def main():
    if not len(sys.argv) == 4:
        msg = """Usage: python get_encode_url.py <Accession Number> <Biological Replicate> <Read>\n
Accession Number\tMust be the one for the entire EXPERIMENT, not for the library, and not for the single fastq file.
Biological Replicate\tInteger number of biological replicate (1, 2, ...)
Read\t\t\tInteger number of single- / paired-end read (1 or 2)\n
Example: python get_encode_url.py ENCSR255XTC 1 2
"""
        print(msg)
        sys.exit(0)
    _, exp_acc, biol_rep, read = sys.argv

    enc_url=get_exp(exp_acc, int(biol_rep), int(read))
    print(enc_url)
    return
     
if __name__ == '__main__':
    main()
