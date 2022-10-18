## Analysis Notebooks

Each Jupyter notebook corresponds to a different aspect of the paper:

- `fitness_assay_simulations.ipynb`

Here, I run simulations of a fitness assays and explore how measurement errors vary with a range of experimental parameters. I also compare these results to approximate theoretical error bounds (obtained assuming Poisson noise in counts)

- `reanalysis_experimental_data.ipynb`

Here, I re-run analysis of the TnSeq data for E. coli B REL606 (the LTEE ancestor) with downsampling and using fewer data points to see how varying sampling/amount of data used for analysis impacts how good measurements are. 

NOTE: `functions.py` contains all the relevant functions for simulating fitness assays and for reanalysis of the sequencing counts data. I import all functions from this script into the notebooks (in order to improve readibility, and not have massive code blocks)


