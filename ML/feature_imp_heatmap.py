#!/usr/bin/env python3
# coding: utf-8

'''
A script to make heatmap of feature importances
'''

struc_features =  ['ncontacts', 'nresidues', 'mech_intra']
struc_features += ['phi_psi', 'sec', 'burr', 'acc']
struc_features += ['IUPRED']

with open('columns_to_consider.txt', 'r', encoding='utf-8') as f:
    for line in f:
        feature = line.rstrip().lstrip()
        if feature in struc_features:
            print(feature)