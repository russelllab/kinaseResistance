# Introduction

This directory contains files and scripts for the analysis of Activark predictions on somatic variants (from [COSMIC](https://cancer.sanger.ac.uk/cosmic/classic) v98).


## Data Files

### COSMIC data

- `data/Cosmic_v98_counts.tsv.gz` </br>
List of all confirmed somatic missense mutations in COSMIC v98, mapped to UniProt protein sequences.

- `data/Cosmic_v98_mismatches.tsv.gz` </br>
List of proteins and number of mismatches when trying to map to Uniprot. These proteins are ignored in the analysis.

These two files are generated using the [cosmic_tools repository](https://github.com/JCGonzS/cosmic_tools).

- `data/Cosmic_v98_cancer_gene_census.csv.gz`</br>
List of genes in the COSMIC Cancer Gene Census. Obtained directly from their data portal.

### Predictor data

- `ML/cosmic_ml_input.txt`</br>
List of all COSMIC kinase somatic variants, used as input for the predictor.

- `ML/cosmic_activark.txt.gz`</br>
Complete predictor output for COSMIC kinase variants.

- `ML/cosmic_activark.log`</br>
Log of predictor output.

#### Analysis results

- `results/ML_output_cosmic_all.tsv.gz`</br>
Final list of kinase somatic variants in COSMIC containing known functional annotations (if available), the predictor scores, the sample counts and the corresponding known roles in cancer for the kinase (if known). Variants are sorted by sample count in descending order.

- `results/ML_output_cosmic_gws.tsv.gz`</br>
Same as above but using only sample counts from genome-wide screens.

- `results/stats_from_cosmic_activark.txt`</br>
Diverse calculations made for `results/ML_output_cosmic_all.tsv.gz`, such us the number of variants predicted as activating, deactivating and/or resistance in the whole dataset or withing known oncogenes or TSGs.

## Scripts

- `merge_cosmic_and_ML.py`</br>
Merges the predictions for somatic variants with COSMIC annotations. Produces `results/ML_output_cosmic_all.tsv.gz` and `results/ML_output_cosmic_gws.tsv.gz`.

- `calculate_stats.py`</br>
Loads results in dataframe and makes a bunch of calculations. 