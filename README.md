# Designing fitness assays

## About

Code corresponding to the following paper:

Resolving deleterious and near-neutral effects requires different pooled fitness assay designs

Anurag Limdi and Michael Baym

In this project, we performed transposon sequencing of the ancestors and evolved clones after 50,000 generations of evolution to identify how the distribution of fitness effects and gene essentiality changes over evolution.

## Organization

### 1. Processed data

This folder is empty: please download the data from https://doi.org/10.5281/zenodo.6547536. This processed data is generated using the scripts in https://github.com/baymlab/2022_Limdi-TnSeq-LTEE (Part 1: Data to Trajectories), and is required for final figure generation and analysis.

### 2. Metadata

This folder contains the relevant metadata for analysis, including gene names, locations, reference genomes, etc.

### 3. Analysis notebooks

Contains the following notebooks:

- fitness_assay_simulations: In this notebook, I estimate theoretical error bounds, and simulate fitness assays using several simplifying assumptions, and explore how modifying experimental parameters impacts errors in measurements.
- reanalysis_experimental_data: In this notebook, I reanalyse our previously published and very deeply sequenced E. coli TnSeq dataset, by downsampling and using a subset of timepoints to explore how changing those parameters impacts measurement errors, and compare results against theory/simulations