#!/usr/bin/env python3

import numpy as np
import scipy as sp
import os, sys, gzip
from sklearn.cluster import KMeans
import seaborn as sns
import pandas as pd
# import tensorflow as tf
# from tensorflow import keras
# from tensorflow.keras import layers
import matplotlib.pyplot as plt
from sklearn import decomposition
from sklearn.preprocessing import MinMaxScaler
import plotly.express as px
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, RepeatedStratifiedKFold
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import cross_val_score
from sklearn.metrics import confusion_matrix, matthews_corrcoef, f1_score, precision_score, recall_score
from sklearn.metrics import auc
from sklearn.metrics import RocCurveDisplay
from sklearn import tree
import xgboost as xgb
import pickle
import argparse

## Old
# 5 3 3 100 LvNA
# 5 3 3 100 AvNL
# 5 4 4 100 AIvNLD
# 5 5 5 50 LDvNAI
# 5 3 5 50 AIvLD
# 5 3 4 100 AvL
# 5 4 4 100 RvN
## New
# 5 5 3 100 RvN
# 5 10 7 100 LDvNAI
# 5 5 7 100 AIvNLD
# 5 7 5 100 AIvLD
# 5 5 3 100 AvL
# 5 10 3 100 AvNL
# 5 12 4 100 LvNA

RANDOM_STATE = 100
ALGO = 'RF' #LR, XGB, RF
N_SPLITS = 10
N_REPEATS = 10
N_JOBS = -1

AA = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

def makeSets(positives, negatives):
    dic = {}
    for set_type, set_type_name in [[positives, 'positives'], [negatives,'negatives']]:
        if set_type_name not in dic:
            dic[set_type_name] = []
        for char in set_type:
            if char == 'A':
                dic[set_type_name].append('activating')
            elif char == 'I':
                dic[set_type_name].append('increase')
            elif char == 'L':
                dic[set_type_name].append('loss')
            elif char == 'D':
                dic[set_type_name].append('decrease')
            elif char == 'R':
                dic[set_type_name].append('resistance')
            elif char == 'N':
                dic[set_type_name].append('neutral')
            else:
                print('Error: invalid set type', char)
                sys.exit(1)
    return dic


def main(max_depth, min_samples_split, min_samples_leaf, n_estimators,\
        name,
        scaler_filename=None, model_filename=None, column_filename=None):
    positives = name.split('v')[0]
    negatives = name.split('v')[1]
    df = pd.read_csv('trainDataFromHitsSplitTrimmedAln.tsv.gz', sep = '\t')
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
                        # 'ATPcount',
                        # 'ChargesWT',
                        # 'ChargesMUT',
                        'ChargesDiff',
                        # 'ATPcount',
                        #   'A_known',
                        #   'D_known',
                        #   'R_known',
                        #   'Phosphomimic',
                        #   'Acetylmimic',
                        #   'hmmScoreWT',
                        #   'hmmScoreMUT',
                          'hmmScoreDiff'
                        ]

    # columns_to_exclude += ['ncontacts', 'nresidues', 'mech_intra']
    # columns_to_exclude += ['phi_psi', 'sec', 'burr', 'acc']
    # columns_to_exclude += ['IUPRED']
    '''
    for aa in AA:
        # if aa not in ['S', 'T', 'Y']:
            columns_to_exclude.append(aa+'_WT')
        # if aa not in ['D', 'E']:
            columns_to_exclude.append(aa+'_MUT')
    '''

    ############
    pfam_ptm_cols = ['ac_pfam', 'me_pfam', 'gl_pfam', 'm1_pfam', 'm2_pfam', 'm3_pfam', 'sm_pfam', 'ub_pfam']
    for i in range(-5,6):
        if i in [-2, -1, 0, 1, 2]: continue
        for col in pfam_ptm_cols:
            columns_to_exclude.append(col.split('_')[0]+'_'+str(i)+'_'+col.split('_')[1])

    pfam_ptm_cols = ['p_pfam']
    for i in range(-5,6):
        if i in [-2, -1, 0, 1, 2]: continue
        for col in pfam_ptm_cols:
            columns_to_exclude.append(col.split('_')[0]+'_'+str(i)+'_'+col.split('_')[1])
    ############

    ptm_cols = ['ac', 'me', 'gl', 'm1', 'm2', 'm3', 'sm', 'ub']
    for i in range(-5,6):
        if i in [-2, -1, 0, 1, 2]: continue
        for col in ptm_cols:
            columns_to_exclude.append(col.split('_')[0]+'_'+str(i))

    ptm_cols = ['p']
    for i in range(-5,6):
        if i in [-2, -1, 0, 1, 2]: continue
        for col in ptm_cols:
            columns_to_exclude.append(col.split('_')[0]+'_'+str(i))

    ############

    adr_cols = ['A', 'D', 'R']
    # adr_cols = ['D', 'R']
    for i in range(-5, 6):
        # if i in [-2, -1, 1, 2]: continue
        for col in adr_cols:
            columns_to_exclude.append(col+'_'+str(i))

    ############

    adr_cols = ['A_pfam', 'D_pfam', 'R_pfam']
    # adr_cols = ['D_pfam', 'R_pfam']
    for i in range(-5, 6):
        # if i in [-2, -1, 0, 1, 2]: continue
        for col in adr_cols:
            columns_to_exclude.append(col.split('_')[0]+'_'+str(i)+'_'+col.split('_')[1])

    df = df.loc[:, ~df.columns.isin(columns_to_exclude)]
    # print (df)

    # scaler = MinMaxScaler()
    # columns_to_scale = ['p_pfam', 'ac_pfam', 'me_pfam', 'gl_pfam', 'm1_pfam', 'm2_pfam', 'm3_pfam', 'sm_pfam', 'ub_pfam']
    # columns_to_scale += ['hmmScoreDiff', 'hmmScoreWT', 'hmmScoreMUT']
    # df[columns_to_scale] = scaler.fit_transform(df[columns_to_scale])

    # print ('columns to consider', df.columns.to_numpy())
    columns_to_consider = '\n'.join(df.columns.to_numpy())
    # print (columns_to_consider)
    if column_filename is not None:
        open(column_filename, 'w').write(columns_to_consider)

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
    dic = makeSets(positives, negatives)
    for row in df.to_numpy():
        # print (row)
        # if row[-1] in ['activating', 'increase']:
        if row[-1] in dic['positives']:
            y.append(1)
            y_names.append(row[-1])
            X.append(row[3:-1])
            train_names.append('/'.join(row[:3]))
        # elif row[-1] in ['loss', 'decrease']:
        elif row[-1] in dic['negatives']:
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
    # print (X)

    X = np.array(X)
    scaler = MinMaxScaler()
    scaler.fit(X)
    X = scaler.transform(X)
    X_test = scaler.transform(X_test)
    # for name, row in zip(test_names, X_test):
    #     if 'A84T' in name:
    #         print (row)
    #         print (len(row))
    #         break
    # sys.exit()
    # pickle.dump(scaler, open('finalized_scaler_RN.pkl', 'wb'))
    if scaler_filename is not None:
        pickle.dump(scaler, open('scaler_'+scaler_filename+'.pkl', 'wb'))

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

    if ALGO == 'LR':
        '''
        parameters = {'C': [0.001, 0.01, 0.1, 1.0, 10.0, 100.0, 1000.0],
                    'solver': ['newton-cg', 'lbfgs', 'liblinear', 'sag', 'saga'],
                    'penalty': ['l1', 'l2', 'elasticnet', 'none'],
                    'max_iter': [100, 250, 500, 1000, 1500, 2000]
                    }
        '''
        parameters = {
                    'C': [0.001, 0.01, 0.1, 1.0, 10.0],
                    'solver': ['lbfgs'],
                    'penalty': ['l2'],
                    'max_iter': [1000]
                    }
        model = LogisticRegression(class_weight='balanced')
        model = GridSearchCV(model, parameters, cv=rskf, scoring='roc_auc', n_jobs=N_JOBS)
        model.fit(X, y)
        clf = LogisticRegression(
                            class_weight='balanced',
                            max_iter=model.best_params_['max_iter'],
                            solver=model.best_params_['solver'],
                            C=model.best_params_['C'],
                            penalty=model.best_params_['penalty']
                            )
    elif ALGO == 'XGB':
        parameters = {
                    'max_depth': [3],
                    'objective': ["binary:logistic"],
                    # 'learning_rate': [0.1],
                    # 'min_samples_split': [2],
                    # 'min_samples_leaf': [3],
                    # 'max_features': ['sqrt', 'log2'],
                    # 'n_estimators': [100]
                    }
        xgb_model = xgb.XGBClassifier(
                        random_state=RANDOM_STATE,
                        scale_pos_weight=float(np.count_nonzero(y == 1))/np.count_nonzero(y == 0)
                        )
        model = GridSearchCV(xgb_model, parameters, cv=rskf, scoring='roc_auc', n_jobs=N_JOBS)
        model.fit(X, y)
        clf = xgb.XGBClassifier(
                # n_estimators=model.best_params_['n_estimators'],
                # leaarning_rate=model.best_params_['learning_rate'],
                # min_samples_leaf=model.best_params_['min_samples_leaf'],
                # min_samples_split=model.best_params_['min_samples_split'],
                max_depth=model.best_params_['max_depth'],
                objective=model.best_params_['objective'],
                # max_features=model.best_params_['max_features'],
                random_state=RANDOM_STATE,
                scale_pos_weight=float(np.count_nonzero(y == 1))/np.count_nonzero(y == 0)
                )
    elif ALGO == 'RF':
        parameters = {
                    'max_depth': max_depth,
                    # 'max_depth': [None],
                    'min_samples_split': min_samples_split,
                    'min_samples_leaf': min_samples_leaf,
                    'criterion': ['gini'],
                    # 'max_features': ['sqrt', 'log2'],
                    'max_features': ['log2'],
                    # 'max_features': ['sqrt'],
                    'n_estimators': n_estimators
                    }
        # parameters = {
        #             'max_depth': [5],
        #             'min_samples_split': [3],
        #             'min_samples_leaf': [3],
        #             'max_features': ['sqrt', 'log2'],
        #             'n_estimators': [100]
        #             }
        # parameters = {
        #             'max_depth': [2,4,5,7,10],
        #             'min_samples_split': [2,3,4,5],
        #             'min_samples_leaf': [2,4,5],
        #             'max_features': ['sqrt', 'log2'],
        #             'n_estimators': [100, 200, 500]
        #             }
        rf = RandomForestClassifier(random_state=RANDOM_STATE, class_weight="balanced", n_jobs=N_JOBS)
        model = GridSearchCV(rf, parameters, cv=rskf, scoring='roc_auc', n_jobs=N_JOBS)
        model.fit(X, y)
        print (name, 'mean_test_score', model.cv_results_['mean_test_score'])
        # print (model.cv_results_['mean_train_score'])
        clf = RandomForestClassifier(
                n_estimators=model.best_params_['n_estimators'],
                min_samples_leaf=model.best_params_['min_samples_leaf'],
                min_samples_split=model.best_params_['min_samples_split'],
                max_depth=model.best_params_['max_depth'],
                max_features=model.best_params_['max_features'],
                criterion=model.best_params_['criterion'],
                random_state=RANDOM_STATE, class_weight="balanced", n_jobs=N_JOBS
                )

    breakLine = '#'.join(['-' for i in range(0, 50)])
    print (breakLine)
    ## Best model hyper-parameters
    print ('Best model found during the CV')
    print (model.best_params_)
    # print (model.predict_proba(X))
    for y_pred, y_true in zip(model.predict_proba(X), y):
        open('ai_ld_roc.txt', 'a').write(str(y_pred[1]) + '\t' + str(y_true) + '\n')
    # sys.exit()
    '''
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
            if ALGO == 'LR':
                clf = LogisticRegression(
                            class_weight='balanced',
                            max_iter=model.best_params_['max_iter'],
                            solver=model.best_params_['solver'],
                            C=model.best_params_['C'],
                            penalty=model.best_params_['penalty']
                            )
            elif ALGO == 'XGB':
                clf = xgb.XGBClassifier(
                    # n_estimators=model.best_params_['n_estimators'],
                    # learning_rate=model.best_params_['learning_rate'],
                    # min_samples_leaf=model.best_params_['min_samples_leaf'],
                    max_depth=model.best_params_['max_depth'],
                    objective=model.best_params_['objective'],
                    # min_samples_split=model.best_params_['min_samples_split'],
                    # max_features=model.best_params_['max_features'],
                    random_state=RANDOM_STATE,
                    scale_pos_weight=float(np.count_nonzero(y == 1))/np.count_nonzero(y == 0)
                    )
            elif ALGO == 'RF':
                clf = RandomForestClassifier(
                    n_estimators=model.best_params_['n_estimators'],
                    min_samples_leaf=model.best_params_['min_samples_leaf'],
                    max_depth=model.best_params_['max_depth'],
                    min_samples_split=model.best_params_['min_samples_split'],
                    max_features=model.best_params_['max_features'],
                    random_state=RANDOM_STATE, class_weight="balanced", n_jobs=N_JOBS
                    )
            clf.fit(X_train, y_train)
            tn, fp, fn, tp = confusion_matrix(y_train, clf.predict(X_train)).ravel()
            #print (tn, fp, fn, tp)
            auc_itr.append(roc_auc_score(y_validation, clf.predict_proba(X_validation)[:,1]))
            mcc_itr.append(matthews_corrcoef(y_validation, clf.predict(X_validation)))
            f1_itr.append(f1_score(y_validation, clf.predict(X_validation)))
            pre_itr.append(precision_score(y_validation, clf.predict(X_validation)))
            rec_itr.append(recall_score(y_validation, clf.predict(X_validation)))
            spe_itr.append(recall_score(y_validation, clf.predict(X_validation), pos_label=0))
            if i == 5:
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
        # print (np.mean(auc_itr), auc_itr)
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
    '''
    ## Fit the best model on the data
    if ALGO == 'LR':
        clf = LogisticRegression(
                    class_weight='balanced',
                    max_iter=model.best_params_['max_iter'],
                    solver=model.best_params_['solver'],
                    C=model.best_params_['C'],
                    penalty=model.best_params_['penalty']
                    )
    elif ALGO == 'XGB':
        clf = xgb.XGBClassifier(
                    # n_estimators=model.best_params_['n_estimators'],
                    # learning_rate=model.best_params_['learning_rate'],
                    # min_samples_leaf=model.best_params_['min_samples_leaf'],
                    max_depth=model.best_params_['max_depth'],
                    objective=model.best_params_['objective'],
                    # max_features=model.best_params_['max_features'],
                    random_state=RANDOM_STATE,
                    scale_pos_weight=float(np.count_nonzero(y == 1))/np.count_nonzero(y == 0)
                    )
    elif ALGO == 'RF':
        clf = RandomForestClassifier(
                    n_estimators=model.best_params_['n_estimators'],
                    min_samples_leaf=model.best_params_['min_samples_leaf'],
                    min_samples_split=model.best_params_['min_samples_split'],
                    max_depth=model.best_params_['max_depth'],
                    max_features=model.best_params_['max_features'],
                    criterion=model.best_params_['criterion'],
                    random_state=RANDOM_STATE, class_weight="balanced", n_jobs=N_JOBS
                    )
    clf.fit(X,y)
    
    if ALGO == 'RF':
        # print (clf.estimator_.decision_path(X))
        '''
        estimator = clf.estimator_
        estimator.fit(X, y)
        text_representation = tree.export_text(estimator)
        # print(text_representation)
        
        fig = plt.figure(figsize=(25,20))
        _ = tree.plot_tree(estimator, 
                        feature_names = feature_names,
                        class_names = y_names,
                        filled=True)
        # plt.show()
        '''
        

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
        df_feature_importances = df_feature_importances.sort_values(by=['Importance'], ascending=False)
        print (df_feature_importances)
        '''
        sns.set(font_scale = 0.6)
        sns.barplot(data=df_feature_importances, color="grey", x="Importance", y="Feature")
        plt.grid(True, lw=0.1)
        '''
        # plt.savefig('feature_imp.png')
        # plt.show()
    

    # filename = 'finalized_model_RN.sav'
    if model_filename is not None:
        pickle.dump(clf, open('model_'+model_filename+'.sav', 'wb'))

    test_types = ['activatingresistance', 'increaseresistance','resistance', 'A', 'TBD', 'Inconclusive', 'TBDincreaseresistance']
    for test_type in test_types:
        print (''.join(['#' for i in range(1,25)]))
        if test_type in ['activatingresistance']:
        # if test_type in ['activatingresistance', 'increaseresistance']:
        # if test_type in ['activatingresistance', 'increaseresistance', 'resistance']:
            X_sub_test = []; y_sub_test = []
            for test_name, p, q in zip(test_names, X_test, y_test):
                if q != test_type: continue
                X_sub_test.append(p)
                y_sub_test.append(1)
            X_sub_test = np.array(X_sub_test)
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
                X_sub_test = np.array(X_sub_test)
                if 'A84' in test_name:
                    print (X_sub_test)
                y_pred = round((clf.predict_proba(X_sub_test)[0])[1], 3)
                print (test_name, y_pred, q)

if __name__ == '__main__':
    '''max_depth = int(sys.argv[1])
    min_samples_split = int(sys.argv[2])
    min_samples_leaf = int(sys.argv[3])
    n_estimators = int(sys.argv[4])'''
    # set arguments
    parser = argparse.ArgumentParser(description='Training for Activark', epilog='End of help.')
    parser.add_argument('max_depth', help='')
    parser.add_argument('min_samples_split', help='')
    parser.add_argument('min_samples_leaf', help='')
    parser.add_argument('n_estimators', help='')
    parser.add_argument('name', help='AIvLD or AIvNLD or LDvNAI or RvN')
    parser.add_argument('--s', help='filename of the scaler to be saved')
    parser.add_argument('--m', help='filename of the model to be saved')
    parser.add_argument('--c', help='filename where the columns to consider must be saved')
    args = parser.parse_args()

    # set input file to default if not provided
    max_depth = int(args.max_depth)
    min_samples_split = int(args.min_samples_split)
    min_samples_leaf = int(args.min_samples_leaf)
    n_estimators = int(args.n_estimators)
    positives = args.name.split('v')[0]
    negatives = args.name.split('v')[1]
    name = args.name

    if args.s: scaler_filename = args.s
    else: scaler_filename = None
    
    if args.m: model_filename = args.m
    else: model_filename = None

    if args.c: column_filename = args.c
    else: column_filename = None
    # print ('hello')
    main([max_depth],[min_samples_split],[min_samples_leaf], [n_estimators],\
        name = name,
        scaler_filename=scaler_filename, model_filename=model_filename, column_filename=column_filename)