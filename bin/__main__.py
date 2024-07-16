import os
import pandas as pd
import numpy as np
import argparse
import time
import import_combine as imp
import perturbation as pert
import perturbation_preprocessing as ppr
from utils import create_logger, get_basename, check_output_dir
from pathlib import Path

# parser.add_argument('--log', action='store_true', help='Enable logging to file')
LOG = True
logger = create_logger(logger_name=__name__, log_file=LOG)

def get_args():
    parser = argparse.ArgumentParser(prog="scPAGE")
    
    parser.add_argument("RawInputPath", default= 'datasets/CellRanger/',
                        help = "File path for Cell Ranger outputs")
    parser.add_argument("Gtf", default='datasets/gtf/hg38.ensGene.gtf',
                        help ="File path for Gtf file")
    parser.add_argument("RiboUrl", default= "http://software.broadinstitute.org/gsea/msigdb/download_geneset.jsp?geneSetName=KEGG_RIBOSOME&fileType=txt",
                        help = "This must be a link to a txt file containing a list of ribosomal RNA")
    parser.add_argument("-n", dest="name", default=None, 
                        help="Output file name without the extension")
    parser.add_argument("-o", dest="output_dir", default="Outputs/", 
                        help="Output directory path")
    parser.add_argument("--subset", dest="sub", action="store_true", 
                        help="If true, the count matrix is reduced in size. Used for testing")
    parser.add_argument("--subset-n", dest="sub-n", default=1000, type=int, 
                        help="Integer for number of cells kept from count matrix. Used for testing")
    parser.add_argument("--subset-sg", dest="sub_sg", default=100, type=int, 
                        help="If true, sgRNA table is reduced")
    parser.add_argument("--umi-threshold", dest="umi_t", default=3, type=int, 
                        help="Minimum threshold for sgRNA UMI")
    parser.add_argument("--subset-sg-n", dest="sub_sg_n", default=1000, type=int, 
                        help="Number of sgRNA to keep")
    parser.add_argument("--convert-ids", dest="conv_ids", action="store_true", 
                        help="If true, convert IDs to gene names")
    parser.add_argument("--id-lookup", dest="id_lookup", default="datasets/Lookup.csv", 
                        help="List of genes for subsetting the sgRNA table")
    parser.add_argument("--percentile-choice", dest="perc_choice", default=0.98, 
                        help="Percentile for counts to keep")
    parser.add_argument("--mito-filter", dest="mito_filter", default=30, 
                        help="Percent threshold for mitochondrial RNA removal")
    parser.add_argument("--ribo-filter", dest="ribo_filter", default=35, 
                        help="Percent threshold for ribosomal RNA removal")

    args = parser.parse_args()

    raw_input = args.RawInputPaht
    gtf = args.gtf
    ribo_url = args.RiboUrl
    coord = args.coordinates
    name = args.name
    output = args.output_dir
    subset = args.subset
    subset_n = args.sub_n
    subset_sg = args.sub_sg
    umi_t = args.umi_t
    subset_sg_n = args.sub_sg_n
    conv_ids = args.conv_ids
    id_l = args.id_lookup
    perc_c = args.perc_choice
    m_filter = args.mito_filter
    r_filter = args.ribo_filter

    return (raw_input, gtf, ribo_url, name, output, subset, subset_n, subset_sg, umi_t, subset_sg_n, conv_ids, id_l, perc_c, m_filter, r_filter)

#RAW_INPUT = path to raw files
if __name__ == "__main__":
    RAW_INPUT, GTF, RIBO_URL, NAME, OUTPUT, SUBSET, SUBSET_N, SUBSET_SG, UMI_T, SUBSET_SG_N, CONV_ID, ID_L, PERCENTILE_CHOICE, MITO_FILTER, RIBO_FILTER = get_args()
    print("Starting Processing Files...")
    time.sleep(2)

    ################### SORT DIRECTORIES ###################
    if NAME is None:
        NAME = get_basename(RAW_INPUT)

    # checks if output exists or create one
    output_dir = check_output_dir(OUTPUT)
    # sets location for plot and table outputs, checks if exists or creates one
    output_plot_dir = check_output_dir(output_dir + '/Plot_outputs')
    output_scanpy_dir = check_output_dir(output_dir + '/Scanpy')
    output_sgRNA_dir = check_output_dir(output_dir + '/sgRNAs')
    
    ################### COMBINE MATRICES ###################
    # set path to raw files
    raw_path = Path(RAW_INPUT)
    parent_raw = str(raw_path.parent.absolute())
    combined_output = check_output_dir(parent_raw + '/Combined/')
    # Combine raw data in scanpy suitable format if not executed yet
    if os.path.exists(combined_output) and os.path.isdir(combined_output):
        if len(os.listdir(combined_output)) > 0:
            print(os.listdir(combined_output))
            print('Combined/ folder is not empty. Please delete old files first')
        else:
            print('Combined/ empty')
            print('Combining Input Matrices')
            # retrieve unique file names
            file_names = set(os.listdir(RAW_INPUT))
            # get unique basenames
            basenames = set(['_'.join(names.split('_')[:-1]) for names in file_names])

            for basename in basenames:
                print(basename)
                imp.combine_inputs(RAW_INPUT, basename, SUBSET, SUBSET_N)
    else:
        print(f'Error. {combined_output} folder does not exist')

    ################### PREPARE SGRNA TABLES ###################
    #  map and subset sgRNAs 
    # check directories exist or create one
    preproc_sgrna = check_output_dir(parent_raw + '/Preprocessing/sgRNAs')
    sg_mapping = check_output_dir(parent_raw + 'sgRNAmapping/')

    for file in os.listdir(sg_mapping):
        print(f'Processing {file}')
        imp.assign_sgRNA_coord(sg_mapping + file, SUBSET_SG, UMI_T, SUBSET_SG_N, preproc_sgrna)

    # get gene IDs for targets and controls
    for files in os.listdir(preproc_sgrna):
        if files.startswith('Coordinates_Targets'):
            imp.sgRNA_geneIDs(preproc_sgrna + files, GTF, ID_L, output_sgRNA_dir, CONV_ID)
        elif files.startswith('Coordinates_Controls'):
            imp.process_ctrl(preproc_sgrna + files, output_sgRNA_dir, CONV_ID)
        else:
            print(f'File {files} not processable')

    ################### PREPROCESSING ###################
    # create output folder
    out_dir_table = check_output_dir(output_scanpy_dir + 'Tables/')
    out_dir_plot = check_output_dir(output_scanpy_dir + 'QC_Plots/')

    for file in os.listdir(combined_output):
        print(file)
        pert.preprocessing(combined_output + file, RIBO_URL, out_dir_plot, out_dir_table, PERCENTILE_CHOICE, MITO_FILTER, RIBO_FILTER)

    ################### PERTURBATION PREPROCESSING ###################
    # list processed sgRNAs
    sgrna_file_names = os.listdir(output_sgRNA_dir)
    # list processed count matrices
    matrix_file_names = os.listdir(out_dir_table)

    # get unique sample identifiers to process controls and targets at the same time
    sample_ids = set(['_'.join(names.split('_')[3:6]) for names in sgrna_file_names])
    matrix_ids = set(['_'.join(names.split('_')[3:6]) for names in matrix_file_names])

    # create list for storing the filtered tables
    all_samples = []
    for id in sample_ids:
        # extract basename
        basename = '_'.join(id.split('_')[1:3])
        # extract sample id
        sample_id = id.split('_')[0]
        # read sgRNA files
        # target
        target_file = [f for f in sgrna_file_names if 'Targets_' + sample_id in f]
        # control
        ctrl_file = [f for f in sgrna_file_names if 'Controls_' + sample_id in f]

        # read equivalent processed matrix
        matrix_file = [f for f in matrix_file_names if basename in f]   
        matrix = pd.read_csv(os.path.join(out_dir_table, matrix_file[0]), index_col=0)

        if len(target_file) == 1 and len(ctrl_file) == 1:
            # read Target sgRNA list
            target = pd.read_csv(os.path.join(output_sgRNA_dir, target_file[0]), sep = '\t')
            # read Control sgRNA list
            ctrl = pd.read_csv(os.path.join(output_sgRNA_dir, ctrl_file[0]), sep = '\t')
            filtered_matrix = ppr.subset_count_table(target, ctrl, matrix, basename, 'datasets/sgRNA_Subset.csv')
            all_samples.append(filtered_matrix)
        else:
            print(f'Error. {target_file}/{ctrl_file} are duplicated files')

        # combine the count tables
        final_combined = ppr.combine_count_table(all_samples)

    ################### PERTURBATION ANALYSIS ###################
    # subset the matrix on PCA
    df_guide, va_genes = pert.most_variable_genes(final_combined,
                        FILTER=True,
                        PERCENTILE = 99)
    # run perturbation significance analysis and generate the graph
    pert.significant_perturbation(df_guide, va_genes, out_dir_table)


