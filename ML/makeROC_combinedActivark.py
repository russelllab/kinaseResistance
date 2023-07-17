#!/usr/bin/env python3
# coding: utf-8

'''
A script to cluster kinase resistance, activating,
and deactivating mutations
'''
import os, sys, gzip
from turtle import position
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.metrics import RocCurveDisplay

def make_curve(name, y_true, y_target, ax):
    viz = RocCurveDisplay.from_predictions(
                                    y_true,
                                    y_target,
                                    name=f"{name}",
                                    alpha=0.75,
                                    lw=2,
                                    ax=ax,
                                )
    # interp_tpr = np.interp(mean_fpr, viz.fpr, viz.tpr)
    # interp_tpr[0] = 0.0
    # tprs.append(interp_tpr)
    # aucs.append(viz.roc_auc)

fig, ax = plt.subplots(figsize=(6, 6))
counter = 1
count_all_map = ['D', 'A']
for files in [
            # 'ai_nld_roc.txt',
            #   'ld_nai_roc.txt',
                'all_roc.txt',
                'all_roc.txt',
                'AIvLD_roc.txt',
                # 'a_nl_roc.txt',
                # 'l_na_roc.txt',
                # 'a_l_roc.txt',
                'RvN_roc.txt',
                ]:
    if 'all_roc' in files:
        counter += 1
    y_true = []
    y_target = []
    for line in open(files):
        name = files.split('_roc')[0].upper()
        if '/' in name: name = name.split('/')[-1]

        if 'all_roc' in files:
            name += '_'+str(count_all_map[counter-2])
            if str(line.split()[-1]) not in ['0','1','2']: continue
            y_target.append(float(line.split()[counter]))
            if float(line.split()[-1].strip())+1 == counter:
                y_true.append(1)
            else:
                y_true.append(0)
            # print (name)
        else:
            if str(line.split()[-1]) not in ['0', '1']: continue
            y_target.append(float(line.split()[2]))
            y_true.append(int(line.split()[3].strip()))
    print (name)
    make_curve(name, y_true, y_target, ax)
    
plt.grid(True)
plt.savefig('roc_all_activark.png', dpi=1000)
plt.savefig('roc_all_activark.svg', dpi=1000)
# plt.show()