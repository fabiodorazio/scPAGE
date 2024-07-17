
**scPAGE**


## Introduction

The purpose of this exercise is to perform a general analysis of the CROP-seq dataset provided by [`Tian et al., 2019`](https://www.sciencedirect.com/science/article/pii/S0896627319306403). In this publication the authors per1formed scRNA-seq combined with CRISPR perturbation during iPSC to neuron transition to interrogate the genes necessary for neuronal survival and functioning.

scPAGE (single cell Perturbation Analysis Gene Expression) aims to perform the necessary steps to generate a list of significantly differentially regulated genes in CRISPR perturbation experiments.
Due to memory and time limitations, adjustments to the pipelines have been made to trim the datasets. It is therefore important to consider this as an exercise and not a comprhensive analysis for research discoveries.
The program generates summary tables and plots in the Output directory

## Quick Start

1. Install [`Conda`](https://conda.io/miniconda.html)

2. Install the required packages from the provided environment.yml file with the following command:

```bash
conda env create -f environment.yml
```

3. cd into the bin directory

4. Start scPAGE with the following command:

```bash
python3 __main.py__
```

## Pipeline Summary

The following are the arguments required for scPAGE:

1. "RawInputPath"           >       "File path for Cell Ranger outputs"
2. "Gtf"                    >       "File path for Gtf file"
3. "RiboUrl"                >       "This must be a link to a txt file containing a list of ribosomal RNA"
4. "-n"                     >       "Output file name without the extension"
5. "-o"                     >       "Output directory path"
6. "--subset"               >       "If true, the count matrix is reduced in size. Used for testing"
7. "--subset-n"             >       "Integer for number of cells kept from count matrix. Used for testing"
8. "--subset-sg"            >       "If true, sgRNA table is reduced"
9. "--umi-threshold"        >       "Minimum threshold for sgRNA UMI"
10. "--subset-sg-n"         >       "Number of sgRNA to keep"
11. "--convert-ids"         >       "If true, convert IDs to gene names"
12. "--id-lookup"           >       "List of genes for subsetting the sgRNA table"
13. "--percentile-choice"   >       "Percentile for counts to keep"
14. "--mito-filter"         >       "Percent threshold for mitochondrial RNA removal"
15. "--ribo-filter"         >       "Percent threshold for ribosomal RNA removal"

## Description and Steps

scPAGE starts by importing the required files provided by the user. It performs some qc checks on the input files:

1. scRNA-seq processing using Scanpy
> Removing of doublets
> Removing cells with low read counts
> Removing lowly represented genes
> Removing cells with too high mitochondrial or ribosomal RNA. (Not stringent) It was difficult to set fixed thresholds as neurons and iPSC had highly variable content of mito and ribo RNAs
> Plot QC metrics
> Generate cleaned count matrix
> Clustering could be performed but it was not necessary for this analysis

2. Perturbation analysis preprocessing
> Select only cells with sgRNA barcode counts > threshold (=100)
> Map sgRNA barcodes to the ref genome GRhg38 sending requests to ['Blat'](https://genome.ucsc.edu/cgi-bin/hgBlat) using the ['requests'](https://pypi.org/project/requests/) module to retrieve ENSEMBL IDs of sgRNAs targets
> Identify control sgRNAs by selecting the non-targeting ones. This is done in absence of clear understanding of control design
> Use the gene names tables to map the ENSEMBL IDs to gene names
> Assign arbitrary gene names to controls (e.i CONTROL1)
> Subset sgRNA table for genes of interest (e.i. rejuvination factors)
> Subset each count table for all samples by cell IDs remaining in the sgRNA table
> Create final table for perturbation analysis:
- Keep only genes expressed in all samples
- Combine conditions and replicates in a single table

3. Perturbation analysis
> Reads are normalised by total number of reads per cell
> The mean read count for cells in the same condition and sgRNA groups is calculated
> Dimensionality reduction is performed with PCA to identify those genes with the highest variance between conditions
> For every condition:
- Fit a normal distribution on the control gene
- Use the parameters of the normal distribution (mean and sd) to test significance between Target and control sgRNA
- Correct p-values using False Discovery Rate
- Plot genes by adjusted p-values
