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

fig, ax = plt.subplots(figsize=(6, 6))
for files in ['ai_nld_roc.txt',
              'ld_nai_roc.txt',
              'ai_ld_roc.txt',
                'a_nl_roc.txt',
                'l_na_roc.txt',
                'a_l_roc.txt',
                'r_n_roc.txt',
                '../data/polyphen2_roc.txt',
                '../data/pmut_roc.txt'
                ]:
    y_target = []
    y_true = []
    for line in open(files):
        y_target.append(float(line.split()[0]))
        y_true.append(float(line.split()[1].strip()))
    name = files.split('_roc')[0].upper()
    if '/' in name:
        name = name.split('/')[-1]
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
plt.grid(True)
plt.savefig('roc_all.png', dpi=1000)
plt.savefig('roc_all.svg', dpi=1000)
# plt.show()