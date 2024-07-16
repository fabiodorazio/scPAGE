import gzip
import logging
import requests
import pandas as pd
import numpy as np

def query_blat(sequence):
    # URL for UCSC BLAT
    url = "https://genome.ucsc.edu/cgi-bin/hgBlat"

    # Payload for the BLAT request
    payload = {
        'type': 'blat',
        'userSeq': sequence,
        'output': 'json'
    }

    # Send the POST request to UCSC BLAT
    response = requests.post(url, data=payload)

    return(response)


def assign_sgRNA_coord(SGRNA,
                SUBSET_N, 
                UMI_T):
    '''
    Formats the sgRNA and filters for UMI count
    Retrieve coordinates of sgRNA targets through Blat
    Assigns coordinates to gene IDs referencing the hg38 gtf file
    '''
    # load the sgRNA file
    sgrna = pd.read_csv(SGRNA, sep = '\t')
    sgrna_file = SGRNA.split('/')[-1].split('.')[0]

    # remove the rows containing the string unprocessed
    sgrna = sgrna[~sgrna['barcode'].str.contains('unprocessed')]
    # remove rows with UMI count less than 3
    sgrna = sgrna[~(sgrna['umi_count'] < UMI_T)]

    # keep only unique sgrnas
    unique_seq = sgrna['barcode'][sgrna['barcode'].duplicated()]

    # query sgRNA sequences through Blast to obtain the target coordinates
    results = []
    for sequence in unique_seq.iloc[:SUBSET_N]:
        print(sequence)
        # query blat
        try:
            r = query_blat(sequence)
            if r.status_code == 200:
                r = r.json()
                # retrieve coordinates for unique and non empty matches
                if r and len(r['blat']) == 1:
                    chr = r['blat'][0][13]
                    start = r['blat'][0][15]
                    end = r['blat'][0][16]
                else:
                    start = 'NA'
                    end = 'NA'

            # retrieve results
            results.append({'barcode':sequence,
                    'chrom':chr,
                    'start':start,
                    'end':end})
        except:
            pass

    # combine final dataframe
    final_df = pd.DataFrame.from_dict(results)
    # remove barcodes with no matches
    final_df.dropna(inplace=True)
    
    # merge with original
    final_df = pd.merge(final_df, sgrna, on='barcode')
    # save
    final_df.to_csv('datasets/Coordinates/' + sgrna_file + '_coordinates.txt', sep = '\t', index = False)

    return(final_df)

import os
path = 'datasets/sgRNAmapping/'
for file in os.listdir(path):
    print(f'processing {file}')
    assign_sgRNA_coord(path + file, 1000, 3)