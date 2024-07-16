import gzip
import logging
import requests
import pandas as pd
import numpy as np
from scipy.io import mmread
from scipy.sparse import csc_matrix
import bioframe as bf


def combine_inputs(PATH, 
                   BASENAME, 
                   SUBSET,
                   N_BEST_COLUMNS):
    '''
    Subsets the input files contained in folder dataset
    Reads in count matrix, cell barcodes and gene names
    Combines the input files in a single table
    '''
    # load count matrix
    logging.info(f'Loading count matrix for {BASENAME}')
    print(PATH + BASENAME + '_matrix.mtx')
    count_matrix = mmread(PATH + BASENAME + '_matrix.mtx')
    # subset the matrix 
    if SUBSET:
        # sum all the values in each column and convert the matrix into a 1d array
        column_sums = np.array(count_matrix.sum(axis=0)).flatten()
        # sort by sum and select the top best columns selected by the user
        col_index = np.argsort(column_sums)[-N_BEST_COLUMNS:]
        # convert the sparse matrix into csr to allow subsetting
        count_matrix = count_matrix.tocsr()
        # subset the matrix
        count_matrix = count_matrix[:,col_index]
        # convert it to array and dataframe
        count_matrix = count_matrix.toarray()
        count_matrix = pd.DataFrame(count_matrix)

    # load genes
    logging.info(f'Loading genes for {BASENAME}')
    genes = pd.read_csv(PATH + BASENAME + '_genes.tsv', sep ='\t', header=None)

    # load cell barcodes
    logging.info(f'Loading cell barcodes for {BASENAME}')
    barcodes = pd.read_csv(PATH + BASENAME + '_barcodes.tsv', sep ='\t', header=None)
    # subset barcodes to match matrix
    barcodes = barcodes.iloc[col_index,:]
    # convert barcodes to list of strings
    barcodes = barcodes[0].astype(str).tolist()

    # combine the datasets
    count_matrix.columns = barcodes
    count_matrix.index = genes[1]

    # subset
    #s_matrix = count_matrix.sample(n=1000, axis=1, random_state=1990)

    # save
    count_matrix.to_csv(PATH + '../Combined/' + BASENAME + '_combined.csv')
    print(count_matrix)

    return(count_matrix)


#combine_inputs('datasets/CellRanger/', 'GSM3543618_iPSC_1', True)


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
                UMI_T,
                READS_T,
                OUT_PATH):
    '''
    Formats the sgRNA and filters for UMI count
    Retrieve coordinates of sgRNA targets through Blat
    Assigns coordinates to gene IDs referencing the hg38 gtf file
    It wasn't clear which were the ctrl sgRNAs therefore I assumed the non-targeting ones as ctrls
    '''
    # load the sgRNA file
    sgrna = pd.read_csv(SGRNA, sep = '\t')
    sgrna_file = SGRNA.split('/')[-1]

    # remove the rows containing the string unprocessed
    sgrna = sgrna[~sgrna['barcode'].str.contains('unprocessed')]
    # remove rows with UMI count < 3
    sgrna = sgrna[~(sgrna['umi_count'] < UMI_T)]
    # remove rows with read counts < 100
    sgrna = sgrna[~(sgrna['read_count'] < READS_T)]

    # query sgRNA sequences through Blast to obtain the target coordinates
    results = []
    for index, row in sgrna.head(SUBSET_N).iterrows():
        cell = row['cell']
        sequence = row['barcode']
        read_count = row['read_count']
        umi_count = row['umi_count']

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
                    start = np.nan
                    end = np.nan

            # retrieve results
            results.append({'cell':cell,
                    'barcode':sequence,
                    'read_count':read_count,
                    'umi_count':umi_count,
                    'chrom':chr,
                    'start':start,
                    'end':end})
        except:
            pass

    # combine final dataframe
    final_df = pd.DataFrame.from_dict(results)
    # save the non-targeting (control) sgRNAs in a separate table
    sg_controls = final_df[final_df.isna().any(axis=1)]

    # remove control barcodes with no matches and save. This is to speed up the gene id retrieval
    final_df.dropna(inplace=True)
    
    # save
    final_df.to_csv(OUT_PATH + 'Coordinates_Targets_' + sgrna_file, sep = '\t', index = False)
    sg_controls.to_csv(OUT_PATH + 'Coordinates_Controls_' + sgrna_file, sep = '\t', index = False)

    return(final_df)
'''
import os
sg_path = 'datasets/sgRNAmapping/'
for file in os.listdir(sg_path):
    print(file)
    assign_sgRNA_coord(sg_path + file, 1000, 3, 50)
'''

def sgRNA_geneIDs(SGRNA,
                  GTF,
                  GENES,
                  OUT_PATH,
                  CONVERT_IDS=True):
    '''
    Retrieves the sequence coordinates of each sgRNA in the dataset
    Matches the coordinates with the gene IDs in the gtf file provided
    Creates an output format compatible with downstream preprocessing
    Retrieves the ensembl Ids and gene names
    Maps out the ensembl Ids in sgRNA to gene names
    '''
    # get file name
    sgrna_file = SGRNA.split('/')[-1]
    SGRNA = pd.read_csv(SGRNA, sep = '\t')

    # read gtf file and rename the columns
    gtf = pd.read_csv(GTF, sep = '\t')
    gtf_columns = ['chrom',
                   'source',
                   'feature',
                   'start',
                   'end',
                   'score',
                   'strand', 
                   'frame',
                   'attribute']
    gtf.columns = gtf_columns

    # ensure to convert coordinates to integers
    SGRNA['start'] = SGRNA['start'].astype(int)
    SGRNA['end'] = SGRNA['end'].astype(int)

    # overlap: bioframe needs the following 3 colnames {chrom, start, end}
    all_overlaps = bf.overlap(SGRNA, gtf)
    
    # extract attribute
    all_overlaps['gene'] = all_overlaps['attribute_'].apply(lambda x:x.split('"')[-2])
                                                            
    # drop columns
    cols_to_keep = ['cell',
                    'barcode',
                    'gene',
                    'read_count']
    all_overlaps = all_overlaps[cols_to_keep]
    # duplicates are generated due to alignment to different parts of the gene. Remove duplicates
    all_overlaps.drop_duplicates(inplace=True)

    # Combine reads from different cells. if duplicated barcodes still exist, they belong to different cells

    if CONVERT_IDS:
        # read gene table
        genes = pd.read_csv(GENES, sep ='\t', header=None)
        # rename columns
        genes.columns = ['gene', 'gene_name']

        # merge
        sgrna_merged = pd.merge(all_overlaps, genes, on='gene')
        # save merged
        sgrna_merged.to_csv(OUT_PATH + 'Filtered_' + sgrna_file, sep = '\t', index = False)
    else:
        # save original
        all_overlaps.to_csv(OUT_PATH + 'Filtered_' + sgrna_file, sep = '\t', index = False)


def process_ctrl(SGRNA,
                 OUT_PATH,
                 CONVERT_IDS=True):
    '''
    Formats the sgRNA ctrl tables
    '''
    # read the file
    sgrna_file = SGRNA.split('/')[-1]
    SGRNA = pd.read_csv(SGRNA, sep = '\t')

    # add gene column
    SGRNA['gene'] = 'Control'
    SGRNA['gene_name'] = 'Control'

    # subset columns
    if CONVERT_IDS:
        cols_to_keep = ['cell',
                        'barcode',
                        'gene',
                        'read_count',
                        'gene_name']
    else:
        cols_to_keep = ['cell',
                        'barcode',
                        'gene',
                        'read_count']
    
    SGRNA = SGRNA[cols_to_keep]
    
    # save
    SGRNA.to_csv(OUT_PATH + 'Filtered_' + sgrna_file, sep = '\t', index = False)


import os
sg_path = 'datasets/Preprocessing/sgRNAs/'
for files in os.listdir(sg_path):
    print(files)
    if files.startswith('Coordinates_Targets'):
        sgRNA_geneIDs(sg_path + files, 'datasets/gtf/hg38.ensGene.gtf', 'datasets/Lookup_IDs.tsv')
    elif files.startswith('Coordinates_Controls'):
        process_ctrl(sg_path + files)
    else:
        print(f'File {files} not processable')





