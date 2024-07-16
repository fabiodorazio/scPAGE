
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
python3 main.py ../assets/ ../assets/genome.fa.gz ../assets/amplicon_coordinates.bed -n Mut_Signatures_Combined_Amplicons --mutational-signatures --combine-amplicons --gene-overlaps
```

## Pipeline Summary

The following 3 arguments are mandatory for scPAGE:

1. "input_amp"    >  "Input amplicon sequences in fasta format"
2. "reference"    >  "Input reference genome in fasta format"
3. "coordinates"  >  "A bed file containing the coordinates of each amplicon"

The following 6 arguments are optional:

1. "-n","name"                  >  "Output file name without the format: can be the name of the analysis, individual, date, etc..."
2. "-qc", "qcmode"              >  "If provided, runs the program in qc-mode only and does not perform the alignment"
3. "-o", est="output_dir"       >  "Output directory path"
4. "--combine-amplicons"        >  "If provided, combines amplicons based on coordinates and basenames"
5. "--mutational-signatures"    >  "If provided, estimates the type of mutational signature"
6. "--gene-overlaps"            >  "Provide gene coordinates to identify the ensembl ID of the gene where the mutation is found"

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
> Subset each count table for all samples by cell IDs remaining is sgRNA table
> Create final table for perturbation analysis:
- Aggregate counts for each gene in cells with same target sgRNA
- Combine conditions and replicates in a single table

3. Perturbation analysis
> Reads are normalised by total number of reads per cell
> The mean read count for cells in the same condition and sgRNA groups is calculated
> Dimensionality reduction is performed with PCA to identify those genes with the highest variance between conditions
> Fit a Gaussian on 


PerturbAnalyser aims to perform a simple mutational signature convertion. This is particularly important for cancer samples where the origin is uncertain.

If target gene/exons coordinates are provided, each mutation can be mapped back and labelled with the gene and exons where it belongs.

Finally, the mutation frequency is calculated for substitutions or signatures and graphics are generated