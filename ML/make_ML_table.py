#!/usr/bin/env python3

import numpy as np
import scipy as sp
import os, sys, gzip
import argparse
import ML

def parse_ML_output(inputFile):
    text = ''
    flag_header = 0
    dic_model = {}
    for line in open(inputFile, 'r'):
        if 'v' in line:
            if line.split('v')[0].isupper() and line.split('v')[1].isupper():
                modelName = line.rstrip().lstrip()
                flagAR = 0
                dic_test = {}
                dic_nld = {}
                dic_params = {}
        if line.startswith('Best model:'):
            params = line.split('{')[1].split('}')[0].split(',')
            for param in params:
                param_name = param.split(':')[0].rstrip().lstrip().replace("'", "")
                param_value = param.split(':')[1].rstrip().lstrip().replace("'", "")
                dic_params[param_name] = param_value
        if line.startswith('MET AUC'):
            metrics = line.rstrip().lstrip().split()[1:]
        if line.startswith('AVG '):
            avgs = line.rstrip().lstrip().split()[1:]
        if line.startswith('STD '):
            stds = line.rstrip().lstrip().split(' ')[1:]
        if line.startswith('activatingresistance results '):
            flagAR = 1
        if line.startswith('REC: ') and flagAR == 1:
            testAR_recall = line.split('REC: ')[1].rstrip().lstrip()
            flagAR = 0
        if line.count('/') == 2:
            dic_test[line.split(' ')[0]] = line.split(' ')[1].rstrip().lstrip()
        if line.startswith('REC: neutral ') or line.startswith('REC: decrease ') or line.startswith('REC: loss '):
            dic_nld[line.split(' ')[1]] = line.split(' ')[2].rstrip().lstrip()
        if line.startswith('REC: decrease '):
            # print (modelName, avgs, testAR_recall, dic_test, dic_nld)
            nld_header = ['neutral', 'loss', 'decrease']
            test_header = ['P46734/MAP2K3/A84T', 'P11309/PIM1/S97N',
                        'O96017/CHEK2/K373E', 'P00533/EGFR/T790M',
                        'P46734/MAP2K3/T222M', 'Q02750/MAP2K1/V211D']
            if flag_header == 0:
                text += 'model\t' + '\t'.join(metrics) + '\t' + 'testAR_recall\t' + '\t'.join(nld_header) + '\t' + '\t'.join(test_header) + '\n'
                flag_header = 1
            text += modelName + '\t'
            for avg in avgs:
                text += avg + '\t'
            text += testAR_recall + '\t'
            for nld in nld_header:
                text += dic_nld[nld] + '\t'
            for test in test_header:
                text += dic_test[test] + '\t'
            text += '\n'
            dic_model[modelName] = dic_params
    print (text)
    print (metrics)
    return dic_model

def make_ML_models(dic_model):
    for model in ['AILDRvN', 'AILDvN', 'AIvNLD', 'AIvN', 'LDvNAI', 'LDvN', 'AIvLD', 'RvN']:
    # for model in ['AIvLD']:
        dic_params = dic_model[model]
        max_depth = [int(dic_params['max_depth'])]
        min_samples_split = [int(dic_params['min_samples_split'])]
        min_samples_leaf = [int(dic_params['min_samples_leaf'])]
        n_estimators = [int(dic_params['n_estimators'])]
        ML.main(max_depth, min_samples_split, min_samples_leaf, n_estimators,
        model, scaler_filename=model, model_filename=model,
        column_filename='columns_to_consider.txt')


if __name__ == '__main__':
    # set arguments
    parser = argparse.ArgumentParser(description='Parse ML output', epilog='End of help.')
    parser.add_argument('inputFile', help='')
    # parser.add_argument('--c', help='filename where the columns to consider must be saved')
    args = parser.parse_args()

    # set input file to default if not provided
    inputFile = args.inputFile

    dic_model = parse_ML_output(inputFile)
    # make_ML_models(dic_model)

