![Website](https://img.shields.io/website?down_message=down&style=plastic&up_message=up&url=http%3A%2F%2Factivark.russelllab.org) ![Twitter Follow](https://img.shields.io/twitter/follow/MechismoNews?style=social)

<img src='webApp/static/img/logo_v2.svg' alt='logo' width='500'/>

Activark is a data-driven, ML-based approach to predict the functional consequence of genetic changes in protein kinases. Activark was trained on a curated dataset of activating (i.e. constitutive-activation or increase in kinase activity), deactivating (i.e. loss or decrease in Kinase activity), and drug-resistance protein variants in human kinases and using sequence and structural features.

Briefly, we applied a random forest algorithm to develop 3 contrasting predictors based on seven types of sequence and structural features:
1. Pred (A v D): The first predictor, activating vs deactivating, represents a typical situation when one has what is believed to be a functional variant (e.g. observed many times in a cohort or dataset) and wishes to distinguish these two possibilities.
2. Pred (A vs D vs N): The second, activating, deactivating or neutral, is more reflective of a situation where one does not know if a variant is functional at all and thus one needs to predict neutrals.
3. Pred (R vs N):The third predictor, resistance vs neutral, predicts if a given mutation is resistant or not.

To access the Activark webservice, go to [here](http://activark.russelllab.org)
To know more about Activark, visit [here](http://activark.russelllab.org/about)

---

## How to run Activark locally?
If you wish run Activark locally on your system, follow the steps below:

Create the environment (activark) from the environment.yml file
> conda env create -f environment.yml

Activate the environment
> conda activate activark

Move to the predictor directory (required)
> cd ML/

Read the help section
> ./prepareTestData.py -h

Example of input:
> ./prepareTestData.py sample_mutations.txt

---

## Performance of Activark

10-fold stratified CV results of all the 3 predictors
![ROC](webApp/static/img/Figure_S3_ML_stats.svg)

---

## Contact

Gurdeep Singh: gurdeep.singh[at]bioquant[dot]uni-heidelberg[dot]de <br>
Torsten Schmenger: torsten.schmenger[at]bioquant[dot]uni-heidelberg[dot]de <br>