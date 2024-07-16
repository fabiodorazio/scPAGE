import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc

def plot_doublets(matrix, basename, output_plot_dir):
    sns.displot(matrix[matrix.prediction == 'doublet'], x = 'difference')
    plt.savefig(f'{output_plot_dir}/{basename}_Doublets_Distr.png')
    plt.close()


    
def plot_qc_metrics(matrix, basename, output_plot_dir):
    '''
    QC Used to get rid of outliers
    Significant higher genes or counts means it is an artifact
    High mitochondrial representation means stressed and apoptotic cells. Set a filter between 5 and 20%
    '''
    
    sc.pl.violin(matrix, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_ribo'], 
             jitter=0.4, multi_panel=True)
    
    plt.savefig(f'{output_plot_dir}/{basename}_Pct_Ribo.png')
    plt.close()



def plot_umap(matrix, output_plot_dir):
    '''
    Plot a UMAP using scanpy
    '''
    sc.pl.umap(matrix)
    plt.savefig(f'{output_plot_dir}/Umap.png')
    plt.close()



def plot_diff_perturb(TABLE,
                      OUT_PLOT_DIR):
    '''
    Plot the adjusted p-values for each gene
    Plot the p-values by condition ordered by the difference between condition
    '''
    # set the style and color palette
    sns.set_style(style="whitegrid")
    palette = sns.color_palette("Set2")

    # set size
    plt.figure(figsize=(12, 8))

    # create the strip plot 
    g = sns.stripplot(data=TABLE, x='gene', y='log_q_value', hue='condition', palette=palette, dodge=True, jitter=True, alpha=0.7, edgecolor='gray')

    # add a horizontal line at -log10(0.05) 
    plt.axhline(-np.log10(0.05), linestyle='--', color='red', linewidth=1.5)

    # rotate the x-axis labels
    plt.xticks(rotation=45, ha='right')

    # set axis labels and title
    plt.xlabel("Gene", fontsize=14)
    plt.ylabel("Log Q-value", fontsize=14)
    plt.title("Differential Gene Expression", fontsize=16)

    # add a legend title and move the legend
    plt.legend(title='Condition', bbox_to_anchor=(1.05, 1), loc='upper left')
    # adjust the layout
    plt.tight_layout()

    # save
    plt.savefig(f'{OUT_PLOT_DIR}/Perturb_Genes.png')
    plt.close()
