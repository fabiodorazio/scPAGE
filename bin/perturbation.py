import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

from scipy.stats import mannwhitneyu, norm, combine_pvalues
from statsmodels.sandbox.stats.multicomp import multipletests
from statsmodels.nonparametric.smoothers_lowess import lowess

import graphic_functions as gfc


def most_variable_genes(MATRIX_PATH,
                        FILTER=True,
                        PERCENTILE = 99):
    """
    Retrieves the 25% cells with most read counts 
    Groups cells by condition and perturbed gene, then calculates the mean
    Performs horizontal dimensionality reduction using PCA
    Selects the 99th percentile of the most variable genes 
    """
    
    # already normalised
    df = pd.read_csv(MATRIX_PATH, index_col=0)

    # create an empty list to contain outputs
    output_list = []

    # trim the dataset
    if FILTER:
        # filter out the bottom 75% covered cells
        # select numeric columns
        num_rows = df.apply(pd.to_numeric, errors='coerce')
        # sum reads for each column
        column_sums = num_rows.sum()
        # sort columns based on sum
        sorted_sums = column_sums.sort_values()
        # select cells with most reads
        top_25_percent_columns = sorted_sums.tail(int(df.shape[1] * 0.25)).index
        df2 = df[top_25_percent_columns]

    else:
        df2 = df
    
    # remove rows with NA
    df2 = df2.dropna()
    
    # group cells by condition and gene, then calculate the mean
    df_t = df2.T.reset_index()
    df_group = df_t.groupby(['index', 'gene']).mean().T
    #df_group.index = df2.index

    # fit PCA for mean and variance then transform to project the original data on the PC axes
    fitted = PCA().fit_transform(df_group)
    r = pd.Series(fitted[:, 1], index=df_group.index).sort_values()

    # Select variable genes based on percentile threshold
    va_genes = r[abs(r) > np.percentile(abs(r), PERCENTILE)].index

    # plot heatmaps
    #plot_VAR(df_group, va_genes)

    # handle output
    output_list.append(df2)
    output_list.append(va_genes)

    return(output_list)


def z_score(x):
        return (x - x.mean()) / x.std()


def significant_perturbation(MATRIX,
                             VA_GENES,
                             OUT_DIR):
    """
    Assess whether a gRNA perturbation is significant
    """

    # initialize an empty DataFrame to store results
    results = pd.DataFrame()
    
    # remove NAs
    MATRIX = MATRIX.dropna()
    # fix first row
    first_row = pd.DataFrame([MATRIX.columns], columns=MATRIX.columns, index=[MATRIX.index.name])
    MATRIX = pd.concat([first_row, MATRIX])

    # subset matrix keeping categorical variables
    new_index_subset = VA_GENES.insert(0,['condition','replicate', 'cell', 'grna', 'gene'])
    # subset, drop columns and transpose
    MATRIX = MATRIX.loc[new_index_subset].dropna()
    MATRIX1 = MATRIX.T
    MATRIX1.drop(columns=['replicate', 'cell'], inplace=True)
    # group by gRNA, reduce to median expression
    df_guide = MATRIX1.groupby(["condition", "gene", "grna"]).median().T

    # combine all control gRNA expression by condition and calculated the median
    try:
        for condition in df_guide.columns.levels[0]:
            ctrl = df_guide.loc[
                :, (df_guide.columns.get_level_values('grna').str.contains("CTRL")) & (df_guide.columns.get_level_values('condition') == condition)
            ].median(axis=1)
        
            # iterate over each individual sgRNA
            for grna in df_guide.columns.levels[2]:
                g = df_guide.loc[
                    :, (df_guide.columns.get_level_values('grna') == grna) & (df_guide.columns.get_level_values('condition') == condition)
                ].squeeze()
                # get the number of cells expressing the grna
                n_cells = MATRIX.loc[:,(MATRIX.loc['grna'] == grna) & (MATRIX.loc['condition'] == condition)].shape[1]

                if g.empty or ctrl.empty:
                    continue

        
                # calculate p-values with non-parametric test (Mann-Whitney U p-value)
                s, p = mannwhitneyu(g.astype(np.float128), ctrl.astype(np.float128))
                # fit a normal distribution on control gRNA expression to extract mean and sd
                params = norm.fit(z_score(ctrl))
                # calculate p-value from fitted distribution. Two sided p-value
                es, ep = combine_pvalues(norm.sf(abs(z_score(g)), *params) * 2)
        
                # append results to the DataFrame
                results = results.append(pd.Series(
                    [condition, g.name[1], grna, s, p, n_cells, es, ep],
                    index=["condition", "gene", "grna", "stat", "p_value", "n_cells", "estat", "ep_value"]
                ), ignore_index=True)
    except:
        pass

    # filter out control gRNAs from the results
    results = results[~results["grna"].str.contains("CTRL")]

    # correct p-values using FDR
    qs = results.groupby(["condition"]
    ).apply(lambda x: pd.Series(multipletests(x['p_value'], method="fdr_bh")[1], index=x['grna'])).reset_index().rename(columns={0: "q_value"})
    results = results.merge(qs)
    results["log_q_value"] = -np.log10(results["q_value"])
    results = results.sort_values("log_q_value", ascending=False)
    results.to_csv(os.path.join(OUT_DIR, "perturbation_assessment.csv"), index=False)

    # plot log_q_values
    gfc.plot_diff_perturb(results)


df_guide = pd.read_csv('Outputs/df_group.csv', index_col=0)
v = pd.read_csv('Outputs/va_genes.csv')
va_genes = pd.Index(v.iloc[:,0], name='condition')
significant_perturbation(df_guide, va_genes, 'Outputs/')

