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
from sklearn import tree

RANDOM_STATE = 0
N_SPLITS = 10
N_REPEATS = 10

AA = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

df = pd.read_csv('trainDataFromTrimmedAln.tsv.gz', sep = '\t')
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
columns_to_exclude = [
                    #'Acc',
                    #'Mutation',
                    #'Gene',
                    'Dataset',
                    'hmmPos',
                    'hmmSS',
                    #   'A_known',
                    #   'D_known',
                    #   'R_known',
                    #   'Phosphomimic',
                    #   'hmmScoreWT',
                    #   'hmmScoreMUT',
                    #   'hmmScoreDiff'
                      ]
for aa in AA:
    if aa not in ['S', 'T', 'Y']:
        columns_to_exclude.append(aa+'_WT')
    if aa not in ['D', 'E']:
        columns_to_exclude.append(aa+'_MUT')

pfam_ptm_cols = ['p_pfam', 'ac_pfam', 'me_pfam', 'gl_pfam', 'm1_pfam', 'm2_pfam', 'm3_pfam', 'sm_pfam', 'ub_pfam']
for i in range(-5,6):
    if i in [0]: continue
    for col in pfam_ptm_cols:
        columns_to_exclude.append(col.split('_')[0]+'_'+str(i)+'_'+col.split('_')[1])

ptm_cols = ['p', 'ac', 'me', 'gl', 'm1', 'm2', 'm3', 'sm', 'ub']
for i in range(-5,6):
    if i in [0]: continue
    for col in pfam_ptm_cols:
        columns_to_exclude.append(col.split('_')[0]+'_'+str(i))

adr_cols = ['A', 'D', 'R']
for i in range(-5, 6):
    if i in [0]: continue
    for col in adr_cols:
        columns_to_exclude.append(col+'_'+str(i))

df = df.loc[:, ~df.columns.isin(columns_to_exclude)]

# scaler = MinMaxScaler()
# columns_to_scale = ['p_pfam', 'ac_pfam', 'me_pfam', 'gl_pfam', 'm1_pfam', 'm2_pfam', 'm3_pfam', 'sm_pfam', 'ub_pfam']
# columns_to_scale += ['hmmScoreDiff', 'hmmScoreWT', 'hmmScoreMUT']
# df[columns_to_scale] = scaler.fit_transform(df[columns_to_scale])

print (df.columns.to_numpy())
# sys.exit()
feature_names = df.columns.to_numpy()
feature_names = feature_names[3:-1]
# sys.exit()
X = []
y = []
train_names = []
y_names = []

X_test = []
y_test = []
test_names = []
for row in df.to_numpy():
    if row[-1] == 'A':
        y.append(1)
        y_names.append(row[-1])
        X.append(row[3:-1])
        train_names.append('/'.join(row[:3]))
    elif row[-1] in ['N', 'D']:
        y.append(0)
        y_names.append(row[-1])
        X.append(row[3:-1])
        train_names.append('/'.join(row[:3]))
    # elif row[-1] == 'R':
    #     y.append(2)
    #     X.append(row[:-1])
    # elif row[-1] != 'R':
    # elif row[-1] == 'Inconclusive':
    else:
        y_test.append(row[-1])
        X_test.append(row[3:-1])
        test_names.append('/'.join(row[:3]))
        # print ('/'.join(row[:3]), row[-1])
        # if row[2] == 'K373E':
        #     print (row)
            # sys.exit()

# for (test_name, y_pred) in zip(test_names, y_test):
#     print (test_name, y_pred)
# sys.exit()

X = np.array(X)
X = np.array(X)
scaler = MinMaxScaler()
scaler.fit(X)
X = scaler.transform(X)
X_test = scaler.transform(X_test)

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
skf = StratifiedKFold(n_splits=N_SPLITS, shuffle=True)
rskf = RepeatedStratifiedKFold(n_splits=N_SPLITS, n_repeats=N_REPEATS)

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

parameters = {'max_depth': [None],
            'min_samples_split': [2],
            'min_samples_leaf': [3],
            'max_features': ['sqrt', 'log2'],
            'n_estimators': [25, 50, 100]
            }
model = RandomForestClassifier(random_state=RANDOM_STATE, class_weight="balanced", n_jobs=-1)
# model = LogisticRegression(class_weight='balanced')
model = GridSearchCV(model, parameters, cv=rskf, scoring='roc_auc', n_jobs=-1)
model.fit(X, y)

breakLine = '#'.join(['-' for i in range(0, 50)])
print (breakLine)
## Best model hyper-parameters
print ('Best model found during the CV')
print (model.best_params_)

# clf = LogisticRegression(class_weight='balanced', max_iter=model.best_params_['max_iter'], solver=model.best_params_['solver'], C=model.best_params_['C'], penalty=model.best_params_['penalty'])
clf = RandomForestClassifier(n_estimators=model.best_params_['n_estimators'],
            min_samples_leaf=model.best_params_['min_samples_leaf'],
            min_samples_split=model.best_params_['min_samples_split'],
            max_depth=model.best_params_['max_depth'],
            max_features=model.best_params_['max_features'],
            random_state=RANDOM_STATE, class_weight="balanced", n_jobs=-1
            )

tprs = []
aucs = []
mean_fpr = np.linspace(0, 1, 100)
fig, ax = plt.subplots(figsize=(6, 6))

AUC= []; MCC= []; F1=[]; PRE=[]; REC=[]; SPE=[]
for i in range(0,10):
    skf = StratifiedKFold(n_splits=N_SPLITS, shuffle=True)
    rskf = RepeatedStratifiedKFold(n_splits=N_SPLITS, n_repeats=N_REPEATS)
    auc_itr = []; mcc_itr= []; f1_itr=[]; pre_itr=[]; rec_itr=[]; spe_itr=[]
    for fold, (train_index, test_index) in enumerate(skf.split(X, y)):
        X_train, X_validation = X[train_index], X[test_index]
        y_train, y_validation = y[train_index], y[test_index]
        # clf = LogisticRegression(class_weight='balanced', max_iter=model.best_params_['max_iter'], solver=model.best_params_['solver'], C=model.best_params_['C'], penalty=model.best_params_['penalty'])
        clf = RandomForestClassifier(n_estimators=model.best_params_['n_estimators'],
            min_samples_leaf=model.best_params_['min_samples_leaf'],
            max_depth=model.best_params_['max_depth'],
            min_samples_split=model.best_params_['min_samples_split'],
            max_features=model.best_params_['max_features'],
            random_state=RANDOM_STATE, class_weight="balanced", n_jobs=-1
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
ax.plot([0, 1], [0, 1], "k--", label="random (AUC = 0.5)")
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
            max_features=model.best_params_['max_features'],
            random_state=RANDOM_STATE, class_weight="balanced", n_jobs=-1
            )
clf.fit(X,y)
# print (clf.estimator_.decision_path(X))
estimator = clf.estimator_
estimator.fit(X, y)
text_representation = tree.export_text(estimator)
# print(text_representation)
fig = plt.figure(figsize=(25,20))
_ = tree.plot_tree(estimator, 
                   feature_names = feature_names,
                   class_names = y_names,
                   filled=True)
plt.show()

print (''.join(['#' for i in range(1,25)]))
data = []
for feature_name, importance in zip(feature_names, clf.feature_importances_):
    if importance > 0.01:
        print (feature_name, importance)
        row = []
        row.append(feature_name)
        row.append(importance)
        data.append(row)

df_feature_importances = pd.DataFrame(data, columns=['Feature', 'Importance'])
df_feature_importances.sort_values(by=['Importance'], ascending=False)
sns.set(font_scale = 0.6)
sns.barplot(data=df_feature_importances, color="grey", x="Importance", y="Feature")
plt.grid(True, lw=0.1)
plt.savefig('feature_imp.png')
# plt.show()

test_types = ['AR', 'R', 'Activating', 'TBD', 'Inconclusive']
for test_type in test_types:
    print (''.join(['#' for i in range(1,25)]))
    if test_type in ['AR', 'R']:
        X_sub_test = []; y_sub_test = []
        for test_name, p, q in zip(test_names, X_test, y_test):
            if q != test_type: continue
            X_sub_test.append(p)
            y_sub_test.append(1)
            
        # print (test_name, round(y_pred[1], 2), y_known)
        print (test_type, 'results', '(', len(X_sub_test), ')')
        # print(roc_auc_score(y_sub_test, clf.predict_proba(X_sub_test)[:,1]))
        # print('MCC:', matthews_corrcoef(y_sub_test, clf.predict(X_sub_test)))
        # print('F1:', f1_score(y_sub_test, clf.predict(X_sub_test)))
        # print('PRE:', precision_score(y_sub_test, clf.predict(X_sub_test)))
        print('REC:', round(recall_score(y_sub_test, clf.predict(X_sub_test)), 3))
        # print('SPE:', recall_score(y_sub_test, clf.predict(X_sub_test), pos_label=0))
    else:
        for test_name, p, q in zip(test_names, X_test, y_test):
            if q != test_type: continue
            X_sub_test = []
            X_sub_test.append(p)
            y_pred = round((clf.predict_proba(X_sub_test)[0])[1], 3)
            print (test_name, y_pred, q)