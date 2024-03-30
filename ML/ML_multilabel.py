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
from sklearn.svm import SVC
from sklearn.naive_bayes import GaussianNB
from sklearn.naive_bayes import MultinomialNB
from sklearn.neural_network import MLPClassifier
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

def log_model(model):
    print ('mean_test_score', model.cv_results_['mean_test_score'])
    # mean of the array
    print (round(model.cv_results_['mean_test_score'].mean(), 2))
    mlflow.log_metric('AUC', round(model.cv_results_['mean_test_score'].mean(), 2))
    # std of the array
    print (round(model.cv_results_['mean_test_score'].std(), 3))
    mlflow.log_metric('AUC_STD', round(model.cv_results_['mean_test_score'].std(), 3))
    # print (model.cv_results_['mean_train_score'])

    breakLine = '#'.join(['-' for i in range(0, 50)])
    print (breakLine)
    ## Best model hyper-parameters
    print ('Best model:', model.best_params_)
    # print (model.predict_proba(X))
    mlflow.log_params(model.best_params_)

def equalNumOfSamples(X, y):
    num_1 = np.count_nonzero(y == 1)
    num_0 = np.count_nonzero(y == 0)
    num_2 = np.count_nonzero(y == 2)
    if num_1 > num_0 and num_2 > num_0:
        # select all rows from X that have y equals 0
        X0 = np.array([X[i] for i in range(0, len(X)) if y[i] == 0])
        # randomly select equal rows as y = 1 from X that have y equals 0
        X1 = np.array([X[i] for i in range(0, len(X)) if y[i] == 1])
        arr_X1 = np.random.choice(len(X1), num_0, replace=False)
        X1 = X1[arr_X1]
        X2 = np.array([X[i] for i in range(0, len(X)) if y[i] == 2])
        arr_X2 = np.random.choice(len(X2), num_0, replace=False)
        X2 = X2[arr_X2]
        X_new = np.concatenate((X0, X1, X2), axis=0)
        y_new = np.concatenate((np.zeros(len(X0)), np.ones(len(X1)), 2*np.ones(len(X2))), axis=0)
    elif num_0 > num_1 and num_2 > num_1:
        # select all rows from X that have y equals 1
        X1 = np.array([X[i] for i in range(0, len(X)) if y[i] == 1])
        # randomly select equal rows as y = 0 from X that have y equals 1
        X0 = np.array([X[i] for i in range(0, len(X)) if y[i] == 0])
        arr_X0 = np.random.choice(len(X0), num_1, replace=False)
        X0 = X0[arr_X0]
        # randomly select equal rows as y = 0 from X that have y equals 1
        X2 = np.array([X[i] for i in range(0, len(X)) if y[i] == 2])
        arr_X2 = np.random.choice(len(X2), num_1, replace=False)
        X2 = X2[arr_X2]
        # print (X0)
        # print (X1)
        X_new = np.concatenate((X0, X1, X2), axis=0)
        y_new = np.concatenate((np.zeros(len(X0)), np.ones(len(X1)), 2*np.ones(len(X2))), axis=0)
    else:
        # select all rows from X that have y equals 1
        X2 = np.array([X[i] for i in range(0, len(X)) if y[i] == 2])
        # randomly select equal rows as y = 0 from X that have y equals 1
        X0 = np.array([X[i] for i in range(0, len(X)) if y[i] == 0])
        arr_X0 = np.random.choice(len(X0), num_2, replace=False)
        X0 = X0[arr_X0]
        # randomly select equal rows as y = 0 from X that have y equals 1
        X1 = np.array([X[i] for i in range(0, len(X)) if y[i] == 2])
        arr_X1 = np.random.choice(len(X1), num_2, replace=False)
        X1 = X1[arr_X1]
        # print (X0)
        # print (X1)
        X_new = np.concatenate((X0, X1, X2), axis=0)
        y_new = np.concatenate((np.zeros(len(X0)), np.ones(len(X1)), 2*np.ones(len(X2))), axis=0)    
    return X_new, y_new

# def main(max_depth, min_samples_split, min_samples_leaf, n_estimators,\
#         name,
#         scaler_filename=None, model_filename=None, column_filename=None):
def main(name, algo='RF',
        scaler_filename=None, model_filename=None, column_filename=None, **kwargs):
    # positives = name.split('v')[0]
    # negatives = name.split('v')[1]
    df = pd.read_csv('trainDataFromHitsSplitTrimmedAln.tsv.gz', sep = '\t')
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
    # dic = makeSets(positives, negatives)
    for row in df.to_numpy():
        # print (row)
        # if row[-1] in ['activating', 'increase']:
        if row[-1] in ['neutral', 'loss', 'decrease', 'increase', 'activating']:
                y_test.append(row[-1])
                X_test.append(row[3:-1])
                test_names.append('/'.join(row[:3]))
        if row[-1] in ['increase', 'activating']:
            y.append(2)
            y_names.append(row[-1])
            X.append(row[3:-1])
            train_names.append('/'.join(row[:3]))
        elif row[-1] in ['decrease', 'loss']:
            y.append(1)
            y_names.append(row[-1])
            X.append(row[3:-1])
            train_names.append('/'.join(row[:3]))
        elif row[-1] in ['neutral']:
            y.append(0)
            y_names.append(row[-1])
            X.append(row[3:-1])
            train_names.append('/'.join(row[:3]))
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
    if scaler_filename is not None:
        # if the algo not in the directory, create it
        if not os.path.exists('scalers/'+algo):
            os.makedirs('scalers/'+algo)
        pickle.dump(scaler, open('scalers/'+algo+'/scaler_'+scaler_filename+'.pkl', 'wb'))

    y = np.array(y)

    ## stratified CV
    # skf = StratifiedKFold(n_splits=N_SPLITS, shuffle=True)
    rskf = RepeatedStratifiedKFold(n_splits=N_SPLITS, n_repeats=N_REPEATS)

    ## To perform the randomizationt test (Salzberg test), enable the this line
    if kwargs['Salzberg'] == 'True':
        np.random.shuffle(y)

    if algo == 'GNB':
        parameters = {
                    'var_smoothing': kwargs['var_smoothing'],
                    }
        model = GaussianNB()
        model = GridSearchCV(model, parameters, cv=rskf, scoring='roc_auc_ovr_weighted', n_jobs=N_JOBS)
        # assign sample weight based on the class distribution
        # like the class_weight parameter in sklearn
        sample_weight = np.zeros(len(y))
        sample_weight[y == 0] = 1.0/np.count_nonzero(y == 0)
        sample_weight[y == 1] = 1.0/np.count_nonzero(y == 1)
        sample_weight[y == 2] = 1.0/np.count_nonzero(y == 2)
        model.fit(X, y, sample_weight=sample_weight)
        print (name)
        log_model(model)
        clf = GaussianNB(
                        var_smoothing=model.best_params_['var_smoothing']
                        )
        clf.fit(X,y, sample_weight=sample_weight)
    elif algo == 'MLP':
        parameters = {
                    'hidden_layer_sizes':kwargs['hidden_layer_sizes'],
                    'activation':kwargs['activation'],
                    'solver':kwargs['solver'],
                    'alpha':kwargs['alpha'],
                    'batch_size':kwargs['batch_size'],
                    'learning_rate':kwargs['learning_rate'],
                    'learning_rate_init':kwargs['learning_rate_init'],
                    'max_iter':kwargs['max_iter'],
                    'shuffle':kwargs['shuffle'],
                    'tol':kwargs['tol'],
                    }
        model = MLPClassifier(random_state=RANDOM_STATE)
        model = GridSearchCV(model, parameters, cv=rskf, scoring='roc_auc_ovr_weighted', n_jobs=N_JOBS)
        # find out if num of 1 is more than 0 or not
        X_new, y_new = equalNumOfSamples(X, y)
        # model.fit(X, y)
        model.fit(X_new, y_new)

        print (name)
        log_model(model)

        clf = MLPClassifier(
                        hidden_layer_sizes=model.best_params_['hidden_layer_sizes'],
                        activation=model.best_params_['activation'],
                        solver=model.best_params_['solver'],
                        alpha=model.best_params_['alpha'],
                        batch_size=model.best_params_['batch_size'],
                        learning_rate=model.best_params_['learning_rate'],
                        learning_rate_init=model.best_params_['learning_rate_init'],
                        max_iter=model.best_params_['max_iter'],
                        shuffle=model.best_params_['shuffle'],
                        tol=model.best_params_['tol'],
                        random_state=RANDOM_STATE
                        )
        clf.fit(X_new,y_new)
    elif algo == 'SVC':
        parameters = {
                    'C': kwargs['C'],
                    'kernel': kwargs['kernel'],
                    'max_iter': [1000],
                    }
        model = SVC(class_weight='balanced', probability=True)
        model = GridSearchCV(model, parameters, cv=rskf, scoring='roc_auc_ovr_weighted', n_jobs=N_JOBS)
        model.fit(X, y)
        
        print (name)
        log_model(model)
        for y_pred, y_true in zip(model.predict_proba(X), y):
            open(name+'_roc.txt', 'w').write(str(y_pred[1]) + '\t' + str(y_true) + '\n')
        # sys.exit()
        
        
        clf = SVC(
                class_weight='balanced',
                max_iter=model.best_params_['max_iter'],
                C=model.best_params_['C'],
                probability=True
                )
        # pass
        clf.fit(X,y)
    elif algo == 'XGB':
        parameters = {
                    'max_depth': kwargs['max_depth'],
                    # 'objective': ["binary:logistic"],
                    # 'learning_rate': [0.1],
                    'min_samples_split': kwargs['min_samples_split'],
                    'min_samples_leaf': kwargs['min_samples_leaf'],
                    'criterion': ['friedman_mse'],
                    'max_features': ['log2'],
                    'n_estimators': kwargs['n_estimators']
                    }
        xgb_model = GradientBoostingClassifier(
                        # n_jobs=N_JOBS,
                        random_state=RANDOM_STATE,
                        # scale_pos_weight=float(np.count_nonzero(y == 1))/np.count_nonzero(y == 0)
                        )
        model = GridSearchCV(xgb_model, parameters, cv=rskf, scoring='roc_auc_ovr_weighted', n_jobs=N_JOBS)
        # define sample weight based on the class distribution
        # like the class_weight parameter in sklearn
        sample_weight = np.zeros(len(y))
        sample_weight[y == 0] = 1.0/np.count_nonzero(y == 0)
        sample_weight[y == 1] = 1.0/np.count_nonzero(y == 1)
        sample_weight[y == 2] = 1.0/np.count_nonzero(y == 2)
        # print (y)
        # print (sample_weight)
        # sys.exit()
        model.fit(X, y, sample_weight=sample_weight)

        log_model(model)

        clf = GradientBoostingClassifier(
                n_estimators=model.best_params_['n_estimators'],
                # leaarning_rate=model.best_params_['learning_rate'],
                min_samples_leaf=model.best_params_['min_samples_leaf'],
                min_samples_split=model.best_params_['min_samples_split'],
                max_depth=model.best_params_['max_depth'],
                # objective=model.best_params_['objective'],
                max_features=model.best_params_['max_features'],
                criterion=model.best_params_['criterion'],
                random_state=RANDOM_STATE,
                # scale_pos_weight=float(np.count_nonzero(y == 1))/np.count_nonzero(y == 0),
                # n_jobs=N_JOBS
                )
        
        clf.fit(X,y, sample_weight=sample_weight)
    elif algo == 'RF':
        parameters = {
                    'max_depth': kwargs['max_depth'],
                    # 'max_depth': [None],
                    'min_samples_split': kwargs['min_samples_split'],
                    'min_samples_leaf': kwargs['min_samples_leaf'],
                    'criterion': ['gini'],
                    # 'max_features': ['sqrt', 'log2'],
                    # 'max_features': ['log2'],
                    'max_features': ['sqrt'],
                    'n_estimators': kwargs['n_estimators']
                    }
        rf = RandomForestClassifier(random_state=RANDOM_STATE, class_weight="balanced", n_jobs=N_JOBS)
        model = GridSearchCV(rf, parameters, cv=rskf, scoring='roc_auc_ovr_weighted', n_jobs=N_JOBS)
        # model = rf
        model.fit(X, y)
        print (name)
        log_model(model)
        for y_pred, y_true in zip(model.predict_proba(X), y):
            open(name+'_roc.txt', 'w').write(str(y_pred[1]) + '\t' + str(y_true) + '\n')
        # sys.exit()
        
        
        clf = RandomForestClassifier(
                    n_estimators=model.best_params_['n_estimators'],
                    min_samples_leaf=model.best_params_['min_samples_leaf'],
                    min_samples_split=model.best_params_['min_samples_split'],
                    max_depth=model.best_params_['max_depth'],
                    max_features=model.best_params_['max_features'],
                    criterion=model.best_params_['criterion'],
                    random_state=RANDOM_STATE, class_weight="balanced", n_jobs=N_JOBS
                    )
        # pass
        clf.fit(X,y)
    
    if algo == 'RF':
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
            # print (feature_name, importance)
            row = []
            row.append(feature_name)
            row.append(importance)
            data.append(row)

        df_feature_importances = pd.DataFrame(data, columns=['Feature', 'Importance'])
        df_feature_importances = df_feature_importances.sort_values(by=['Importance'], ascending=False)
        # print (df_feature_importances)
        df_feature_importances.to_csv('feature_imp_'+name+'.csv', index=False)
        
        # sns.set(font_scale = 0.6)
        # sns.barplot(data=df_feature_importances, color="grey", x="Importance", y="Feature")
        # plt.grid(True, lw=0.1)
        
        # plt.savefig('feature_imp.png')
        # plt.savefig('feature_imp_'+name+'.svg', format='svg', dpi=600)
        # sys.exit()
        # plt.show()
    

    # filename = 'finalized_model_RN.sav'
    if model_filename is not None:
        # if the algo not in the directory, create it
        if not os.path.exists('models/'+algo):
            os.makedirs('models/'+algo)
        pickle.dump(clf, open('models/'+algo+'/model_'+model_filename+'.sav', 'wb'))
        # pickle.dump(clf, open('models/model_'+model_filename+'.sav', 'wb'))

    test_types = ['activatingresistance',
                'increaseresistance',
                'resistance',
                'A',
                'TBD',
                'Inconclusive',
                'TBDincreaseresistance',
                'TBDloss',
                'activating',
                'increase',
                'neutral',
                'loss',
                'decrease']
    # test_types = ['activatingresistance', 'increaseresistance','resistance', 'A', 'TBD', 'Inconclusive', 'TBDincreaseresistance']
    for test_type in test_types:
        print (''.join(['#' for i in range(1,25)]))
        if test_type in ['activatingresistance']:
        # if test_type in ['activatingresistance', 'increaseresistance']:
        # if test_type in ['activatingresistance', 'increaseresistance', 'resistance']:
            X_sub_test = []; y_sub_test = []
            names = []
            for test_name, p, q in zip(test_names, X_test, y_test):
                if q != test_type: continue
                X_sub_test.append(p)
                y_sub_test.append(1)
                names.append(test_name)
            X_sub_test = np.array(X_sub_test)
            # print (test_name, round(y_pred[1], 2), y_known)
            print (test_type, 'results', '(', len(X_sub_test), ')')
            # print(roc_auc_score(y_sub_test, clf.predict_proba(X_sub_test)[:,1]))
            # print('MCC:', matthews_corrcoef(y_sub_test, clf.predict(X_sub_test)))
            # print('F1:', f1_score(y_sub_test, clf.predict(X_sub_test)))
            # print('PRE:', precision_score(y_sub_test, clf.predict(X_sub_test)))
            for name, y_prob, y_pred in zip(names, clf.predict(X_sub_test), clf.predict_proba(X_sub_test)):
                print (name, y_prob, y_pred)
            y_pred = [1 if i in [2] else 0 for i in clf.predict(X_sub_test)]
            print('REC:', recall_score(y_sub_test,
                                        # clf.predict(X_sub_test),
                                        y_pred,
                                        # pos_label=['2'],
                                        # labels=['2', '1', '0'],
                                        # average='weighted'
                                        ))
            # print('SPE:', recall_score(y_sub_test, clf.predict(X_sub_test), pos_label=0))
            # store in MLflow
            mlflow.log_metric(test_type+'_REC', round(recall_score(y_sub_test, y_pred), 3))
        else:
            pred_neutral = []; known_neutral = []
            pred_deactivating = []; known_deactivating = []
            pred_activating = []; known_activating = []
            for test_name, p, q in zip(test_names, X_test, y_test):
                if q != test_type: continue
                X_sub_test = []
                X_sub_test.append(p)
                X_sub_test = np.array(X_sub_test)
                # if 'A84' in test_name:
                #     print (X_sub_test)
                # print ((clf.predict_proba(X_sub_test)[0]))
                y_outcome = clf.predict(X_sub_test)[0]
                y_pred = round((clf.predict_proba(X_sub_test)[0])[2], 3)
                y_pred_deact = round((clf.predict_proba(X_sub_test)[0])[1], 3)
                y_pred_neutral = round((clf.predict_proba(X_sub_test)[0])[0], 3)
                if 'L858R' in test_name:
                    print (test_name, clf.predict_proba(X_sub_test)[0], clf.predict(X_sub_test)[0], q)
                if q == 'neutral':
                    known_neutral.append(1)
                    # pred_neutral.append(1 if y_pred_neutral>=0.5 else 0)
                    pred_neutral.append(1 if y_outcome==0 else 0)
                if q in ['loss', 'decrease']:
                    known_deactivating.append(1)
                    # pred_deactivating.append(1 if y_pred_deact>=0.5 else 0)
                    pred_deactivating.append(1 if y_outcome in [1] else 0)
                if q in ['increase', 'activating']:
                    known_activating.append(1)
                    # pred_activating.append(1 if y_pred>=0.5 else 0)
                    pred_activating.append(1 if y_outcome in [2] else 0)
                if q in ['resistance', 'neutral', 'loss', 'decrease', 'increase', 'activating']:
                    continue
                # print (clf.predict_proba(X_sub_test)[0])
                print (test_name, clf.predict_proba(X_sub_test)[0], clf.predict(X_sub_test)[0], q)
            if test_type == 'neutral':
                print('REC:', test_type, round(recall_score(known_neutral, pred_neutral), 3))
                mlflow.log_metric(test_type+'_REC', round(recall_score(known_neutral, pred_neutral), 3))
            if test_type in ['loss', 'decrease']:
                print('REC:', test_type, round(recall_score(known_deactivating, pred_deactivating), 3))
                mlflow.log_metric(test_type+'_REC', round(recall_score(known_deactivating, pred_deactivating), 3))
            if test_type in ['increase', 'activating']:
                print('REC:', test_type, round(recall_score(known_activating, pred_activating), 3))
                mlflow.log_metric(test_type+'_REC', round(recall_score(known_activating, pred_activating), 3))

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
    # positives = args.name.split('v')[0]
    # negatives = args.name.split('v')[1]
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
    print ('END')
