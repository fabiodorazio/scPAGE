import pandas as pd
import scanpy as sc
import scvi
import graphic_functions as gpf
import numpy as np


# read matrix and transpose it to align with scanpy requirements: rows as individual cells
#COMBINED_MATRIX = 'datasets/Combined/GSM3543618_iPSC_1_combined.csv'

#SCDATA = sc.read_csv(COMBINED_MATRIX).T

def filter_out_doublets(COMBINED_MATRIX, 
                        DOUBLETS):
    '''
    Uses the labelled doublet cells to filter them out of the original matrix
    '''
    # reload the raw matrix
    SCDATA = sc.read_csv(COMBINED_MATRIX).T
    # make obs and var unique if gene duplication exists
    SCDATA.var_names_make_unique()
    SCDATA.obs_names_make_unique()

    # add doublet column with true or false values based on the match with the doublets dataset
    SCDATA.obs['doublet'] = SCDATA.obs.index.isin(DOUBLETS.index)
    # filter out predicted doublets from the original table
    singlet_matrix = SCDATA[~SCDATA.obs.doublet]
    
    return(singlet_matrix)


def label_doublets(SCDATA,
                   BASENAME,
                   OUT_PLOT):
    '''
    Preprocesses the data by removing lowly expressed genes
    Selects the top highly variable genes to feed the model
    Removing doublets cells by employing the scvi training algorithm
    Creates an output table that 
    '''
    # Clean the dataset for training the model
    # make obs and var unique if gene duplication exists
    SCDATA.var_names_make_unique()
    SCDATA.obs_names_make_unique()
    # filter only genes with expression in at least 10 cells
    sc.pp.filter_genes(SCDATA, min_cells=10)
    # keep the top highly variable genes that describe the model at best
    sc.pp.highly_variable_genes(SCDATA, n_top_genes=2000, subset=True, flavor='seurat_v3')
    # train model to predict the doublets
    scvi.model.SCVI.setup_anndata(SCDATA) # model setup with default parameters
    v = scvi.model.SCVI(SCDATA)
    v.train()
     # create and train the SOLO model to predict doublets
    solo = scvi.external.SOLO.from_scvi_model(v)
    solo.train()
    # predictions on whether a cell is a doublet or not. The higher the score the more likely the cell is a doublet 
    predictions = solo.predict()
    # make a new column as predicted label
    predictions['prediction'] = solo.predict(soft=False)

    # count doublets and singlets. Expected values are between 5 and 20% of doublets
    singlets = predictions.groupby('prediction').count()
    print(f"number of singles it {singlets.iloc[1,1]}")
    # calculate the difference between doublets and singlets to be less stringent on the filter for doublets
    predictions['difference'] = predictions.doublet - predictions.singlet
    # plot the count distribution on the predicted doublets
    gpf.plot_doublets(predictions, BASENAME, OUT_PLOT)
    # predict doublets for cells that have label doublets and a difference higher than 1
    doublets = predictions[(predictions.prediction == 'doublet') & (predictions.difference > 0.5)]

    return(doublets)



def label_mithocondrial_genes(SCDATA):
    '''
    Mitochondrial genes are labelled using MT- at the start
    '''
    SCDATA.var['mt'] = SCDATA.var.index.str.startswith('MT-')

    return(SCDATA)


def label_ribosomal_genes(SCDATA, RIBO_URL): # OPTIONAL
    '''
    Import list of ribosomal genes from broad institute
    Filter out genes
    '''
    # import list of ribosomal genes in pandas from broad institute
    ribo_genes = pd.read_table(RIBO_URL, skiprows=2, header=None)
    # add column for ribosomal genes on original matrix 
    SCDATA.var['ribo'] = SCDATA.var_names.isin(ribo_genes[0].values)

    return(SCDATA)


def calculate_qc_metrics(SCDATA, 
                         RIBO_URL,
                         BASENAME,
                         OUT_PLOT):
    '''
    Use scanpy to calculate qc metrics
    '''
    # label mitochondrial and ribosomal RNAs
    mito_matrix = label_mithocondrial_genes(SCDATA)
    ribo_mito_matrix = label_ribosomal_genes(mito_matrix, RIBO_URL)

    # the new columns added mt and ribo are fed to the function to calculate the qc
    sc.pp.calculate_qc_metrics(ribo_mito_matrix, qc_vars=['mt', 'ribo'], percent_top=None, log1p=False, inplace=True)

    # plot
    gpf.plot_qc_metrics(ribo_mito_matrix, BASENAME, OUT_PLOT)
    # remove genes that are not in at least 4 cells
    sc.pp.filter_genes(ribo_mito_matrix, min_cells=4)

    # filter by total counts. Arbitrary based on the data
    sc.pp.filter_cells(ribo_mito_matrix, min_genes=200)

    return(ribo_mito_matrix)



def filtering(SCDATA, 
              PERCENTILE_CHOICE, 
              MITO_FILTER, 
              RIBO_FILTER):
    '''
    Filter out cells based on the number of genes by counts
    Mitocondrial and ribosomal rna content
    '''
    # calculate upper limit with the 98th percentile. You can set an arbitrary number of genes. Ex = 3000
    U_limit = np.quantile(SCDATA.obs.n_genes_by_counts.values, PERCENTILE_CHOICE)
    # filter genes by U_limit
    SCDATA = SCDATA[SCDATA.obs.n_genes_by_counts < U_limit]

    # filter out cells that have mitochondrial rna reads > 10%
    SCDATA = SCDATA[SCDATA.obs.pct_counts_mt < MITO_FILTER]
    # filter out cells that have ribosomal rna reads > 35%
    SCDATA = SCDATA[SCDATA.obs.pct_counts_ribo < RIBO_FILTER]
    
    return(SCDATA)


def normalisation(SCDATA):
    '''
    sc sequencing the variation is high between cells
    Normalisation is required for comparison between cells
    There is often an order of magnitude difference between the counts per cell in the dataset
    '''
    # normalize every cell to 10,000 UMI so that the total number fo counts per cell is equal
    sc.pp.normalize_total(SCDATA, target_sum=1e4)
    # change to log counts
    sc.pp.log1p(SCDATA)
    # save the new matrix in the raw table in matrix. This populates the raw table
    SCDATA.raw = SCDATA
    return(SCDATA)


# all together
def preprocessing(COMBINED_MATRIX, 
                  RIBO_URL,
                  OUT_PLOT,
                  OUT_TABLE,
                  PERCENTILE_CHOICE,
                  MITO_FILTER,
                  RIBO_FILTER):
    '''
    Combines functions to generate a filtered count matrix
    '''
    # get basename
    basename = COMBINED_MATRIX.split('/')[-1].split('.')[0]
    # read matrix and transpose it to align with scanpy requirements: rows as individual cells
    input_matrix = sc.read_csv(COMBINED_MATRIX).T
    
    # label and remove doublets
    doublets = label_doublets(input_matrix, basename, OUT_PLOT)
    singlets = filter_out_doublets(COMBINED_MATRIX, doublets)
    
    # calculate qc metrics
    SCDATA = calculate_qc_metrics(singlets, RIBO_URL, basename, OUT_PLOT)
    
    SCDATA = filtering(SCDATA, PERCENTILE_CHOICE, MITO_FILTER, RIBO_FILTER)

    # convert to pd dataframe for saving
    output_df = pd.DataFrame(SCDATA.X, index=SCDATA.obs_names, columns=SCDATA.var_names).T
    # save table
    output_df.to_csv(OUT_TABLE + basename + '_Preprocessed.csv')

    return(SCDATA)
