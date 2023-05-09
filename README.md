![Website](https://img.shields.io/website?down_message=down&style=plastic&up_message=up&url=http%3A%2F%2Factivark.russelllab.org) ![Twitter Follow](https://img.shields.io/twitter/follow/MechismoNews?style=social)

![ActivarkLogo](webApp/static/img/ActivarkLogo.png)

Activark is a predictor of activating, deactivating and resistance mutations in kinases
Access the Flask webApp [here](http://activark.russelllab.org)

---

## How to run Activark?

Move to the directory (mandatory)
> cd ML/

See the help section
> ./prepareTestData.py -h

Example:
> ./prepareTestData.py sample_mutations.txt

---

## Performance

10-fold stratified CV for the AD predictor
![ROC](ML/roc_CV.png)

10-fold stratified CV on shuffled dataset
![ROC_shuffled](ML/roc_CV_shuffled.png)

---

## Feature importance

![Performance](ML/feature_imp.png)

---

## Best Decision tree

![Performance](ML/estimator.png)

---

## Database

Activark uses the following DB schema:
![Database](DB/schema/diagrams/summary/relationships.implied.compact.png)