# Activark predictions on variants from UniProt

## Data Files

### UniProt data

- `data/allUniProtVariants.txt`</br>
List of all UniProt protein variants and their annotations.

### Predictor data

- `ML/uniprot_muts_activark.txt.gz`</br>
Complete predictor output for UniProt variants.

- `ML/cosmic_activark.log`</br>
Log of predictor output.

### Analysis results

- `results/uniprot_muts_activark_with_annotations.txt.gz`</br>
Final list of kinase variants in COSMIC containing known functional annotations (if available), the predictor scores, UniProt annotations, number of Pubmed IDs the variant has been reported in, and the genetic disease acronym (if found).

- `results/stats_from_uniprot_activark.tsv`</br>
Diverse calculations made for different subsets of variants in `results/uniprot_muts_activark_with_annotations.txt.gz`, such us the fractions of variants predicted as activating, deactivating and/or resistance.

- `results/stats_from_uniprot_activark_cols.tsv`</br>
Column descriptions for `results/stats_from_uniprot_activark.tsv`.

- `results/stats_from_uniprot_diseases_all.tsv`</br>
Diverse numbers calculated for all UniProt diseases with kinase variants.

- `results/stats_from_uniprot_diseases_unkvars.tsv`</br>
Same but only for disease that have uncharacterized variants.

- `results/stats_from_uniprot_diseases_cols.tsv`</br>
Column descriptions for `results/stats_from_uniprot_diseases_all.tsv` and `...unkvars.tsv`.

### Plots

Contains diverse scatter plots showing genetic diseases distribution according to the fractions of variants predicted as activating or deactivating. Dots are coloured according to the fuctional outcome of characterized variants.  

## Scripts
- `merge_and_calculate_stats.py`</br>
Merges the predictions for UniProt variants and their annotations. Produces all the 'Analysis results' files.
