import os
import pandas as pd
import numpy as np

def genes_of_interest(TARGET,
                      CTRL,
                      GENES_INTEREST):
    '''
    Subset sgRNA tables on gene of interest
    '''
    # read genes of interest table
    genes_of_interest = pd.read_csv(GENES_INTEREST, sep = ',')

    # keep only relevant genes
    sgrna_subset = TARGET[TARGET['gene_name'].isin(genes_of_interest['gene_name'])]
 
    # append controls
    sgrna_subset_ctrl = pd.concat([sgrna_subset, CTRL], axis=0)

    return(sgrna_subset_ctrl)


def normalise(matrix):
    """
    Normalize expression by total number of counts per cell.
    """
    # divide counts by total number of read counts per cell. Log2 transform
    norm_matrix = np.log2(1 + matrix.apply(lambda x: x / float(x.sum()), axis=0) * 1e4)

    return(norm_matrix)


###
def subset_count_table(TARGET,
                       CTRL,
                       MATRIX,
                       BASENAME,
                       GENES_INTEREST):
    '''
    Uses the sgRNA table to extract the gene of interest
    The sgRNA-associated count matrix is subset on cells carrying the guide for those selected genes  
    Labels groups by condition and replicate

    '''
    # subset genes of interest
    goInt = genes_of_interest(TARGET, CTRL, GENES_INTEREST)

    # subset matrix based on gene of interest
    #cells = goInt['cell'].tolist()
     # only columns in df
    #col_subset = [col for col in cells if col in MATRIX.columns]
    # subset df
    #matrix_subset = MATRIX[col_subset]

    # normalise read counts 
    norm_matrix = normalise(MATRIX)

    # create an empty dataframe to fill
    df_matrix = pd.DataFrame()
    for gene_name in goInt['gene_name'].unique():
        # slice the list of sgRNA
        matrix_slice = goInt[goInt['gene_name'] == gene_name]
        # retrieve cells
        cells_gene = matrix_slice['cell']
        print(type(cells_gene))
        try:
            # keep only columns with cell barcodes that match the sgRNA
            matrix_subset = norm_matrix[cells_gene]
            # create a new row with perturbed gene labels
            new_row_gene = pd.DataFrame([[gene_name] * len(cells_gene)], columns=cells_gene, index=['gene'])
            # create a new row with sample
            new_row_sample = pd.DataFrame([[BASENAME.split('_')[0]] * len(cells_gene)], columns=cells_gene, index=['condition'])
            # create a new row with replicate
            new_row_repl = pd.DataFrame([['rep' + BASENAME.split('_')[1]] * len(cells_gene)], columns=cells_gene, index=['replicate'])

            # concatenate rows with matrix
            new_matrix_subset = pd.concat([new_row_sample, 
                                           new_row_repl, 
                                           new_row_gene, 
                                           matrix_subset])
        except:
            pass
        
        # concatenate all
        df_matrix = pd.concat([df_matrix, new_matrix_subset], axis=1)
        # add index to columns
        col_index = pd.Index(df_matrix.columns, name='cell')
        df_matrix.columns = col_index
        
    return(df_matrix)


### LIST_DF = all_samples
def combine_count_table(LIST_DF):
    '''
    Loads the count matrix
    Keeps only the common genes and sorts them alphabetically
    Loads the sgRNA tables
    Binds Target and Controls tables
    Retrieves only the cell barcodes found in the sgrna table. These are those cells were a gene KO has been detected
    '''
    # list unique genes 
    common_genes = set(LIST_DF[0].index)
    for df in LIST_DF[1:]:
        common_genes.intersection_update(df.index)

    # keep only common genes
    common_genes = list(common_genes)

    # sort the dfs indexes alphabetically
    new_list_df_sorted = [df.loc[common_genes] for df in LIST_DF]

    # bind all tables together
    all_together = pd.concat(new_list_df_sorted, axis=1)

    return(all_together)


