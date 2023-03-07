import numpy as np
import scipy as sp
import os, sys, gzip
from sklearn.cluster import KMeans
import seaborn as sns
import pandas as pd
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
import matplotlib.pyplot as plt
from sklearn import decomposition
from sklearn.preprocessing import MinMaxScaler
import plotly.express as px
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, RepeatedStratifiedKFold
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score
from sklearn.metrics import confusion_matrix, matthews_corrcoef, f1_score, precision_score, recall_score
from sklearn.metrics import auc
from sklearn.metrics import RocCurveDisplay

AA = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

df = pd.read_csv('trainData.tsv.gz', sep = '\t')
df['Dataset'] = df['Dataset'].replace(to_replace='train', value=0.025, regex=True)
df['Dataset'] = df['Dataset'].replace(to_replace='test', value=0.3, regex=True)
# exclude columns
# df = df.loc[:, ~df.columns.isin(['allHomologs','exclParalogs','specParalogs','orthologs', 'bpso','bpsh'])]
df = df.loc[:, ~df.columns.isin([
                            # 'allHomologs',
                            # 'exclParalogs',
                            # 'specParalogs',
                            # 'orthologs'
                            # 'bpso',
                            # 'bpsh'
                            ])]
# exclude columns to make the data matrix
original_df = df.copy()
columns_to_exclude = ['Acc',
                      'Mutation',
                      'Gene',
                      'Dataset',
                      'hmmPos',
                      'hmmSS',
                    #   'A_known',
                    #   'D_known',
                    #   'R_known',
                    #   'Phosphomimic',
                    #   'hmmScoreWT',
                    #   'hmmScoreMUT',
                      'hmmScoreDiff'
                      ]
# for aa in AA:
#     columns_to_exclude.append(aa+'_WT')
#     columns_to_exclude.append(aa+'_MUT')
pfam_ptm_cols = ['p_pfam', 'ac_pfam', 'me_pfam', 'gl_pfam', 'm1_pfam', 'm2_pfam', 'm3_pfam', 'sm_pfam', 'ub_pfam']
for i in range(-5,6):
    for col in pfam_ptm_cols:
        columns_to_exclude.append(col.split('_')[0]+'_'+str(i)+'_'+col.split('_')[1])

df = df.loc[:, ~df.columns.isin(columns_to_exclude)]

# scaler = MinMaxScaler()
# columns_to_scale = ['p_pfam', 'ac_pfam', 'me_pfam', 'gl_pfam', 'm1_pfam', 'm2_pfam', 'm3_pfam', 'sm_pfam', 'ub_pfam']
# columns_to_scale += ['hmmScoreDiff', 'hmmScoreWT', 'hmmScoreMUT']
# df[columns_to_scale] = scaler.fit_transform(df[columns_to_scale])

print (df.columns.to_numpy())
feature_names = df.columns.to_numpy()
feature_names = feature_names[:-1]
# sys.exit()
X = []
y = []

X_test = []
y_test = []
for row in df.to_numpy():
    if row[-1] == 'A':
        y.append(1)
        X.append(row[:-1])
    elif row[-1] == 'D':
        y.append(0)
        X.append(row[:-1])
    # elif row[-1] == 'R':
    #     y.append(2)
    #     X.append(row[:-1])
    # elif row[-1] == 'R':
    #     y_test.append(0)
    #     X_test.append(row[:-1])

X = np.array(X)
X = np.array(X)
scaler = MinMaxScaler()
scaler.fit(X)
X = scaler.transform(X)
# X_test = scaler.transform(X_test)

y = np.array(y)

'''
skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=0)
skf.get_n_splits(X, y)

print(skf)
clf = RandomForestClassifier(max_depth=2, random_state=0, class_weight="balanced")
print (cross_val_score(clf, X, y, scoring="roc_auc", cv=skf))
# clf.fit(X, y)
# for y_pred, y_known in zip(clf.predict(X_test), y_test):
#     print (y_pred, y_known)
sys.exit()

for i, (train_index, test_index) in enumerate(skf.split(X, y)):
    print(f"Fold {i}:")
    print (X[train_index])
    sys.exit()
    print(f"  Train: index={train_index}")
    print(f"  Test:  index={test_index}")
'''
## stratified CV
skf = StratifiedKFold(n_splits=10)
rskf = RepeatedStratifiedKFold(n_splits=10, n_repeats=10)

## To perform the randomizationt test (Salzberg test), enable the this line
# np.random.shuffle(y)

parameters = {'C': [0.001, 0.01, 0.1, 1.0, 10.0, 100.0, 1000.0],
            'solver': ['newton-cg', 'lbfgs', 'liblinear', 'sag', 'saga'],
            'penalty': ['l1', 'l2', 'elasticnet', 'none'],
            'max_iter': [100, 250, 500, 1000, 1500, 2000]
            }

parameters = {'C': [0.001],
            'solver': ['lbfgs'],
            'penalty': ['l2'],
            'max_iter': [100, 250, 500]
            }

parameters = {'max_depth': [2],
            'min_samples_leaf': [5],
            'n_estimators': [10]
            }
model = RandomForestClassifier(random_state=0, class_weight="balanced")
# model = LogisticRegression(class_weight='balanced')
model = GridSearchCV(model, parameters, cv=rskf, scoring='roc_auc', n_jobs=-1)
model.fit(X, y)

breakLine = '-'.join(['-' for i in range(0, 50)])
print (breakLine)
## Best model hyper-parameters
print ('Best model found during the CV')
print (model.best_params_)

# clf = LogisticRegression(class_weight='balanced', max_iter=model.best_params_['max_iter'], solver=model.best_params_['solver'], C=model.best_params_['C'], penalty=model.best_params_['penalty'])
clf = RandomForestClassifier(n_estimators=model.best_params_['n_estimators'],
            min_samples_leaf=model.best_params_['min_samples_leaf'],
            max_depth=model.best_params_['max_depth'],
            random_state=0, class_weight="balanced"
            )

tprs = []
aucs = []
mean_fpr = np.linspace(0, 1, 100)
fig, ax = plt.subplots(figsize=(6, 6))

AUC= []; MCC= []; F1=[]; PRE=[]; REC=[]; SPE=[]
for i in range(0,10):
    skf = StratifiedKFold(n_splits=10, shuffle=True)
    rskf = RepeatedStratifiedKFold(n_splits=5, n_repeats=10)
    auc_itr = []; mcc_itr= []; f1_itr=[]; pre_itr=[]; rec_itr=[]; spe_itr=[]
    for fold, (train_index, test_index) in enumerate(skf.split(X, y)):
        X_train, X_validation = X[train_index], X[test_index]
        y_train, y_validation = y[train_index], y[test_index]
        # clf = LogisticRegression(class_weight='balanced', max_iter=model.best_params_['max_iter'], solver=model.best_params_['solver'], C=model.best_params_['C'], penalty=model.best_params_['penalty'])
        clf = RandomForestClassifier(n_estimators=model.best_params_['n_estimators'],
            min_samples_leaf=model.best_params_['min_samples_leaf'],
            max_depth=model.best_params_['max_depth'],
            random_state=0, class_weight="balanced"
            )
        clf.fit(X_train, y_train)
        tn, fp, fn, tp = confusion_matrix(y_train, model.predict(X_train)).ravel()
        #print (tn, fp, fn, tp)
        auc_itr.append(roc_auc_score(y_validation, model.predict_proba(X_validation)[:,1]))
        mcc_itr.append(matthews_corrcoef(y_validation, model.predict(X_validation)))
        f1_itr.append(f1_score(y_validation, model.predict(X_validation)))
        pre_itr.append(precision_score(y_validation, model.predict(X_validation)))
        rec_itr.append(recall_score(y_validation, model.predict(X_validation)))
        spe_itr.append(recall_score(y_validation, model.predict(X_validation), pos_label=0))
        if i == 0:
            viz = RocCurveDisplay.from_estimator(
                                                clf,
                                                X_validation,
                                                y_validation,
                                                name=f"ROC fold {fold}",
                                                alpha=0.3,
                                                lw=1,
                                                ax=ax,
                                            )
            interp_tpr = np.interp(mean_fpr, viz.fpr, viz.tpr)
            interp_tpr[0] = 0.0
            tprs.append(interp_tpr)
            aucs.append(viz.roc_auc)
    AUC.append(np.mean(auc_itr))
    # print (np.mean(auc_itr))
    MCC.append(np.mean(mcc_itr))
    F1.append(np.mean(f1_itr))
    PRE.append(np.mean(pre_itr))
    SPE.append(np.mean(spe_itr))
    REC.append(np.mean(rec_itr))

#####################################################
ax.plot([0, 1], [0, 1], "k--", label="chance level (AUC = 0.5)")
mean_tpr = np.mean(tprs, axis=0)
mean_tpr[-1] = 1.0
mean_auc = auc(mean_fpr, mean_tpr)
std_auc = np.std(aucs)
ax.plot(
    mean_fpr,
    mean_tpr,
    color="b",
    label=r"Mean ROC (AUC = %0.2f $\pm$ %0.2f)" % (mean_auc, std_auc),
    lw=2,
    alpha=0.8,
)

std_tpr = np.std(tprs, axis=0)
tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
ax.fill_between(
    mean_fpr,
    tprs_lower,
    tprs_upper,
    color="grey",
    alpha=0.2,
    label=r"$\pm$ 1 std. dev.",
)

ax.set(
    xlim=[-0.05, 1.05],
    ylim=[-0.05, 1.05],
    xlabel="False Positive Rate",
    ylabel="True Positive Rate",
    title=f"Mean ROC curve with variability",
)
ax.axis("square")
ax.legend(loc="lower right")
plt.show()
#####################################################
    

print (np.std(AUC))
## Cross-validation results
breakLine = '-'.join(['-' for i in range(0, 50)])
print (breakLine)
print ('Stratified CV results')
print ('MET', 'AUC ', 'MCC ', 'F1  ', 'PRE ', 'REC ', 'SPE')
print ('AVG', round(np.mean(AUC),2),round(np.mean(MCC),2),round(np.mean(F1),2),round(np.mean(PRE),2),round(np.mean(REC),2),round(np.mean(SPE),2))
print ('STD', round(np.std(AUC),3),round(np.std(MCC),3),round(np.std(F1),3),round(np.std(PRE),3),round(np.std(REC),3),round(np.std(SPE),3))
print ('Number of act mutations in the train set:', np.count_nonzero(y))
print ('Number of deact mutations in the train set:', len(y) - np.count_nonzero(y))

## Fit the best model on the data
# clf = LogisticRegression(class_weight='balanced', max_iter=model.best_params_['max_iter'], solver=model.best_params_['solver'], C=model.best_params_['C'], penalty=model.best_params_['penalty'])
clf = RandomForestClassifier(n_estimators=model.best_params_['n_estimators'],
            min_samples_leaf=model.best_params_['min_samples_leaf'],
            max_depth=model.best_params_['max_depth'],
            random_state=0, class_weight="balanced"
            )
clf.fit(X,y)
for feature_name, importance in zip(feature_names, clf.feature_importances_):
    if importance > 0: print (feature_name, importance)