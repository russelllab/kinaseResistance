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
from sklearn.svm import SVC
from sklearn.naive_bayes import GaussianNB
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.ensemble import VotingClassifier
from sklearn.model_selection import cross_val_score
from sklearn.metrics import confusion_matrix, matthews_corrcoef, f1_score, precision_score, recall_score
from sklearn.metrics import auc
from sklearn.metrics import RocCurveDisplay
from sklearn import tree
import pickle
import argparse
import mlflow
import mlflow.sklearn
from mlflow.models import infer_signature

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

def equalNumOfSamples(X, y):
    num_1 = np.count_nonzero(y == 1)
    num_0 = np.count_nonzero(y == 0)
    if num_1 > num_0:
        # select all rows from X that have y equals 0
        X0 = np.array([X[i] for i in range(0, len(X)) if y[i] == 0])
        # randomly select equal rows as y = 1 from X that have y equals 0
        X1 = np.array([X[i] for i in range(0, len(X)) if y[i] == 1])
        arr_X1 = np.random.choice(len(X1), num_0, replace=False)
        X1 = X1[arr_X1]
        X_new = np.concatenate((X0, X1), axis=0)
        y_new = np.concatenate((np.zeros(len(X0)), np.ones(len(X1))), axis=0)
    else:
        # select all rows from X that have y equals 1
        X1 = np.array([X[i] for i in range(0, len(X)) if y[i] == 1])
        # randomly select equal rows as y = 0 from X that have y equals 1
        X0 = np.array([X[i] for i in range(0, len(X)) if y[i] == 0])
        arr_X0 = np.random.choice(len(X0), num_1, replace=False)
        X0 = X0[arr_X0]
        # print (X0)
        # print (X1)
        X_new = np.concatenate((X0, X1), axis=0)
        y_new = np.concatenate((np.zeros(len(X0)), np.ones(len(X1))), axis=0)
    
    return X_new, y_new


# def main(max_depth, min_samples_split, min_samples_leaf, n_estimators,
#         name, algo=None,
#         scaler_filename=None, model_filename=None, column_filename=None, **kwargs):
def main(name,
         scaler_filename=None, model_filename=None, column_filename=None, **kwargs):
    positives = name.split('v')[0]
    negatives = name.split('v')[1]
    df = pd.read_csv('trainDataFromHitsSplitTrimmedAln.tsv.gz', sep = '\t', low_memory=False)
    df['Dataset'] = df['Dataset'].replace(to_replace='train', value=0.025, regex=True)
    df['Dataset'] = df['Dataset'].replace(to_replace='test', value=0.3, regex=True)
    # exclude columns
    # df = df.loc[:, ~df.columns.isin(['allHomologs','exclParalogs','specParalogs','orthologs', 'bpso','bpsh'])]
    df = df.loc[:, ~df.columns.isin([
                                # 'allHomologs',
                                # 'exclParalogs',
                                # 'specParalogs',
                                # 'orthologs',
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
                        # 'ChargesWT',
                        # 'ChargesMUT',
                        # 'ChargesDiff',
                        # 'ATPcount',
                        #   'A_known',
                        #   'D_known',
                        #   'R_known',
                        #   'Phosphomimic',
                        'ReversePhosphomimic',
                        #   'Acetylmimic',
                          'ReverseAcetylmimic',
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
    pfam_ptm_cols = ['me_pfam', 'gl_pfam', 'm1_pfam', 'm2_pfam', 'm3_pfam', 'sm_pfam', 'ub_pfam']
    for i in range(-5,6):
        if i in [-2, -1, 0, 1, 2]: continue
        for col in pfam_ptm_cols:
            columns_to_exclude.append(col.split('_')[0]+'_'+str(i)+'_'+col.split('_')[1])

    pfam_ptm_cols = ['p_pfam', 'ac_pfam']
    for i in range(-5,6):
        if i in [-2, -1, 0, 1, 2]: continue
        for col in pfam_ptm_cols:
            columns_to_exclude.append(col.split('_')[0]+'_'+str(i)+'_'+col.split('_')[1])
    ############

    ptm_cols = ['me', 'gl', 'm1', 'm2', 'm3', 'sm', 'ub']
    for i in range(-5,6):
        if i in [-2, -1, 0, 1, 2]: continue
        for col in ptm_cols:
            columns_to_exclude.append(col.split('_')[0]+'_'+str(i))

    ptm_cols = ['p', 'ac']
    for i in range(-5,6):
        if i in [-2, -1, 0, 1, 2]: continue
        for col in ptm_cols:
            columns_to_exclude.append(col.split('_')[0]+'_'+str(i))

    ############

    adr_cols = ['A', 'D', 'R']
    # adr_cols = ['D', 'R']
    for i in range(-5, 6):
        # if i in [-2, -1, 0, 1, 1]: continue
        for col in adr_cols:
            columns_to_exclude.append(col+'_'+str(i))

    ############

    adr_cols = ['A_pfam', 'D_pfam', 'R_pfam']
    # adr_cols = ['D_pfam', 'R_pfam']
    for i in range(-5, 6):
        # if i in [-2, -1, 0, 1, 2]: continue
        for col in adr_cols:
            columns_to_exclude.append(col.split('_')[0]+'_'+str(i)+'_'+col.split('_')[1])

    adr_cols = ['A', 'D', 'R']
    # adr_cols = ['R']
    AA = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 
        'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    for mut_type in adr_cols:
        for aa in AA:
            col = mut_type+'_'+aa
            for i in range(-5, 6):
                if i in [-2, -1, 0, 1, 2]: continue
                columns_to_exclude.append(col+'_'+str(i)+'_pfam')

    hom_cols = ['allHomologs', 'exclParalogs', 'specParalogs', 'orthologs', 'bpso', 'bpsh']
    # hom_cols = ['allHomologs', 'exclParalogs', 'specParalogs', 'orthologs']
    for hom_type in hom_cols:
        for i in range(-5, 6):
            if i in [0]: continue
            columns_to_exclude.append(hom_type+'_'+str(i))
    
    tax_cols = ['eukaryotes', 'mammals', 'metazoa', 'vertebrates']
    for tax_type in tax_cols:
        for i in range(-5, 6):
            if i in [0]: continue
            columns_to_exclude.append(tax_type+'_'+str(i))


    df = df.loc[:, ~df.columns.isin(columns_to_exclude)]
    # print (df)
    
    # scaler = MinMaxScaler()
    # columns_to_scale = ['p_pfam', 'ac_pfam', 'me_pfam', 'gl_pfam', 'm1_pfam', 'm2_pfam', 'm3_pfam', 'sm_pfam', 'ub_pfam']
    # columns_to_scale += ['hmmScoreDiff', 'hmmScoreWT', 'hmmScoreMUT']
    # df[columns_to_scale] = scaler.fit_transform(df[columns_to_scale])

    # print ('columns to consider', df.columns.to_numpy())
    columns_to_consider = '\n'.join(df.columns.to_numpy())
    # print (columns_to_consider)

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
        if row[-1] in ['neutral', 'loss', 'decrease']:
                y_test.append(row[-1])
                X_test.append(row[3:-1])
                test_names.append('/'.join(row[:3]))
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
        pickle.dump(scaler, open('scalers/scaler_'+scaler_filename+'.pkl', 'wb'))

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
    if kwargs['Salzberg'] == 'True':
        np.random.shuffle(y)

    estimators = []
    for algo in ['XGB', 'RF', 'GNB', 'MLP', 'SVC']:
        # mlflow.set_experiment(algo)
        # print all run info in the experiment
        # print (mlflow.get_experiment_by_name(algo))
        # print all runs in the experiment
        df = mlflow.search_runs(experiment_ids=mlflow.get_experiment_by_name(algo).experiment_id)
        # print header names
        # print (df.columns)
        # print runName
        # print (df['tags.mlflow.runName'])
        # get runID of runName = 'all'
        runID = df[df['tags.mlflow.runName'] == name]['run_id'].to_numpy()[0]
        print (runID)
        # get params of the runID
        params = mlflow.get_run(runID).data.params
        print (params)
        if algo == 'GNB':
            model = GaussianNB(var_smoothing=float(params['var_smoothing']))
            estimators.append(('GNB', model))
        elif algo == 'MLP':
            model = MLPClassifier(hidden_layer_sizes=eval(params['hidden_layer_sizes']),
                                activation=params['activation'],
                                solver=params['solver'],
                                alpha=float(params['alpha']),
                                batch_size=params['batch_size'],
                                learning_rate=params['learning_rate'],
                                learning_rate_init=float(params['learning_rate_init']),
                                max_iter=int(params['max_iter']),
                                shuffle=eval(params['shuffle']),
                                tol=float(params['tol']),
                                random_state=RANDOM_STATE)
            estimators.append(('MLP', model))
        elif algo == 'RF':
            model = RandomForestClassifier(max_depth=int(params['max_depth']),
                                        min_samples_split=int(params['min_samples_split']),
                                        min_samples_leaf=int(params['min_samples_leaf']),
                                        n_estimators=int(params['n_estimators']),
                                        class_weight="balanced",
                                        random_state=RANDOM_STATE)
            estimators.append(('RF', model))
        elif algo == 'XGB':
            model = GradientBoostingClassifier(max_depth=int(params['max_depth']),
                                            min_samples_split=int(params['min_samples_split']),
                                            min_samples_leaf=int(params['min_samples_leaf']),
                                            n_estimators=int(params['n_estimators']),
                                            random_state=RANDOM_STATE)
            estimators.append(('XGB', model))
        elif algo == 'SVC':
            model = SVC(C=float(params['C']),
                        kernel=params['kernel'],
                        probability=True,
                        random_state=RANDOM_STATE)
            estimators.append(('SVC', model))
    print (estimators)
    clf = VotingClassifier(estimators=estimators, voting='soft')
    # sample_weight = np.zeros(len(y))
    # sample_weight[y == 0] = 1.0/np.count_nonzero(y == 0)
    # sample_weight[y == 1] = 1.0/np.count_nonzero(y == 1)
    # clf = clf.fit(X, y, sample_weight=sample_weight)
    X_new, y_new = equalNumOfSamples(X, y)
    clf = clf.fit(X_new, y_new)
    # clf = clf.fit(X, y)
    # sys.exit()

    # filename = 'finalized_model_RN.sav'
    if model_filename is not None:
        pickle.dump(clf, open('models/model_'+model_filename+'.sav', 'wb'))
    
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
            clf = VotingClassifier(estimators=estimators, voting='soft')
            X_train_new, y_train_new = equalNumOfSamples(X_train, y_train)
            clf.fit(X_train_new, y_train_new)
            # tn, fp, fn, tp = confusion_matrix(y_train, clf.predict(X_train)).ravel()
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
        title=f"Mean ROC curve with variability of {N_SPLITS} folds in {name}",
    )
    ax.axis("square")
    ax.legend(loc="lower right")
    plt.grid(linewidth=0.5)
    plt.savefig('ensemble_'+runName+'_roc.svg', format='svg', dpi=1000)

    mlflow.log_metric("AUC", round(np.mean(AUC),2))
    mlflow.log_metric("STD_AUC", round(np.std(AUC),3))
    mlflow.log_metric("MCC", round(np.mean(MCC),2))
    mlflow.log_metric("STD_MCC", round(np.std(MCC),3))
    mlflow.log_metric("F1", round(np.mean(F1),2))
    mlflow.log_metric("STD_F1", round(np.std(F1),3))
    mlflow.log_metric("PRE", round(np.mean(PRE),2))
    mlflow.log_metric("STD_PRE", round(np.std(PRE),3))
    mlflow.log_metric("REC", round(np.mean(REC),2))
    mlflow.log_metric("STD_REC", round(np.std(REC),3))
    mlflow.log_metric("SPE", round(np.mean(SPE),2))
    mlflow.log_metric("STD_SPE", round(np.std(SPE),3))
    mlflow.log_metric("N_POS", np.count_nonzero(y))
    mlflow.log_metric("N_NEG", len(y) - np.count_nonzero(y))
    #####################################################

    test_types = ['activatingresistance', 'increaseresistance','resistance', 'A', 'TBD', 'Inconclusive', 'TBDincreaseresistance', 'activating', 'neutral', 'loss', 'decrease']
    # test_types = ['activatingresistance', 'increaseresistance','resistance', 'A', 'TBD', 'Inconclusive', 'TBDincreaseresistance']
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
            # store in MLflow
            mlflow.log_metric(test_type+'_REC', round(recall_score(y_sub_test, clf.predict(X_sub_test)), 3))
            # print('SPE:', recall_score(y_sub_test, clf.predict(X_sub_test), pos_label=0))
        else:
            pred_neutral = []; known_neutral = []
            pred_deactivating = []; known_deactivating = []
            for test_name, p, q in zip(test_names, X_test, y_test):
                if q != test_type: continue
                X_sub_test = []
                X_sub_test.append(p)
                X_sub_test = np.array(X_sub_test)
                # if 'A84' in test_name:
                #     print (X_sub_test)
                if 'L858R' in test_name:
                    print (test_name, clf.predict_proba(X_sub_test)[0], clf.predict(X_sub_test)[0], q)
                y_pred = round((clf.predict_proba(X_sub_test)[0])[1], 3)
                if q == 'neutral':
                    known_neutral.append(1)
                    if q in dic['negatives']:
                        pred_neutral.append(1 if y_pred<0.5 else 0)
                    else:
                        pred_neutral.append(1 if y_pred>=0.5 else 0)
                    label = 1 if q in dic['negatives'] else 0
                    # print (test_name, y_pred, q, label, negatives)
                if q in ['loss', 'decrease']:
                    known_deactivating.append(1)
                    if q in dic['negatives']:
                        pred_deactivating.append(1 if y_pred<0.5 else 0)
                    else:
                        pred_deactivating.append(1 if y_pred>=0.5 else 0)
                if q in ['resistance', 'neutral', 'loss', 'decrease']:
                    continue
                print (test_name, y_pred, q)
            if test_type == 'neutral':
                print('REC:', test_type, round(recall_score(known_neutral, pred_neutral), 3))
                # store in MLflow
                mlflow.log_metric(test_type+'_REC', round(recall_score(known_neutral, pred_neutral), 3))
            if test_type in ['loss', 'decrease']:
                print('REC:', test_type, round(recall_score(known_deactivating, pred_deactivating), 3))
                # store in MLflow
                mlflow.log_metric(test_type+'_REC', round(recall_score(known_deactivating, pred_deactivating), 3))

if __name__ == '__main__':
    mlflow.set_tracking_uri(uri="http://127.0.0.1:8080")
    mlflow.set_experiment('Ensemble')
    for runName in ['AIvLD', 'RvN']:
        with mlflow.start_run(run_name=runName) as run:
            main(name=runName, Salzberg='False')
