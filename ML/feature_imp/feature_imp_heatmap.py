#!/usr/bin/env python3
# coding: utf-8

'''
A script to make heatmap of feature importances
'''

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

MODEL = 'XGB'
INPUT_FILE = 'feature_imp_AIvLD.csv'
AA = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
        'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
VMAX = 0.055
VMIN = 0.0
LINEWIDTH = 0.05

dic_name2value = {}
with open(MODEL+'/'+INPUT_FILE, 'r', encoding='utf-8') as f:
    for line in f:
        if line.startswith('Feature,Importance'): continue
        feature = line.rstrip().split(',')[0]
        value = float(line.rstrip().split(',')[1])
        dic_name2value[feature] = value


struc_features =  ['ncontacts', 'nresidues', 'mech_intra']
struc_features += ['phi_psi', 'sec', 'burr', 'acc']
struc_features += ['IUPRED', 'ATPcount']

misc_features = ['ChargesWT', 'ChargesMUT', 'ChargesDiff']
misc_features += ['Phosphomimic', 'Acetylmimic']
misc_features += ['hmmScoreWT', 'hmmScoreMUT']

one_hotencoding_features = []
for suffix in ['WT', 'MUT']:
    for aa in AA:
        one_hotencoding_features.append(aa + '_' + suffix)

data = []
window_features = []

#####################################
######## one-hot encoding features ########
#####################################

## one-hot encoding features
for feature in one_hotencoding_features:
    window_features.append(feature)
    data.append(dic_name2value[feature])

df = pd.DataFrame(data, columns=['0'], index=window_features)
sns.set(font_scale = 0.3)
sns.heatmap(df,
            cmap=sns.color_palette("Greys", as_cmap=True),
            xticklabels=True,
            yticklabels=True,
            vmin=VMIN,
            vmax=VMAX,
            linewidths=LINEWIDTH,
            linecolor='black',
            )
# plt.show()
plt.savefig(MODEL+'/'+INPUT_FILE.split('.')[0]+'_ohe.svg', dpi=1000)
data = []
window_features = []

#####################################
######## individual features ########
#####################################

## misc features
for feature in misc_features:
    window_features.append(feature)
    data.append(dic_name2value[feature])

## homology features
homology_features = ['allHomologs', 'exclParalogs', 'specParalogs', 'orthologs', 'bpso', 'bpsh']
for feature in homology_features:
    window_features.append(feature)
    data.append(dic_name2value[feature + '_0'])

## subset features
subset_features = ['metazoa', 'vertebrates', 'mammals', 'eukaryotes']
for feature in subset_features:
    window_features.append(feature)
    data.append(dic_name2value[feature + '_0'])

## structure features
for feature in struc_features:
    window_features.append(feature)
    data.append(dic_name2value[feature])

df = pd.DataFrame(data, columns=['0'], index=window_features)
sns.set(font_scale = 0.3)
sns.heatmap(df,
            cmap=sns.color_palette("Greys", as_cmap=True),
            xticklabels=True,
            yticklabels=True,
            vmin=VMIN,
            vmax=VMAX,
            linewidths=LINEWIDTH,
            linecolor='black',
            )
# plt.show()
plt.savefig(MODEL+'/'+INPUT_FILE.split('.')[0]+'_invidual.svg', dpi=1000)
data = []
window_features = []

#####################################
######## window features ########
#####################################

## PTM Features
ptm_features = ['p', 'ub', 'ac', 'm1', 'm2', 'm3', 'me', 'sm', 'gl']
for feature in ptm_features:
    window_features.append(feature)
    row = []
    for i in range(-2, 3):
        row.append(dic_name2value[feature + '_' + str(i)])
    data.append(row)

for feature in ptm_features:
    window_features.append(feature + '_pfam')
    row = []
    for i in range(-2, 3):
        row.append(dic_name2value[feature + '_' + str(i) + '_pfam'])
    data.append(row)

df = pd.DataFrame(data, columns=['-2', '-1', '0', '1', '2'], index=window_features)
sns.set(font_scale = 0.2)
sns.heatmap(df,
            cmap=sns.color_palette("Greys", as_cmap=True),
            xticklabels=True,
            yticklabels=True,
            vmin=VMIN,
            vmax=VMAX,
            linewidths=LINEWIDTH,
            linecolor='black',
            )
# plt.show()
plt.savefig(MODEL+'/'+INPUT_FILE.split('.')[0]+'_window.svg', dpi=1000)
data = []
window_features = []

## MUT_TYPE features
mutType_features = ['A', 'D', 'R']
for feature in mutType_features:
    for aa in AA:
        window_features.append(feature + '_' + aa)
        row = []
        for i in range(-2, 3):
            row.append(dic_name2value[feature + '_' + aa + '_' + str(i) + '_pfam'])
        data.append(row)

# print (data)
df = pd.DataFrame(data, columns=['-2', '-1', '0', '1', '2'], index=window_features)
sns.set(font_scale = 0.2)
sns.heatmap(df,
            cmap=sns.color_palette("Greys", as_cmap=True),
            xticklabels=True,
            yticklabels=True,
            vmin=VMIN,
            vmax=VMAX,
            linewidths=LINEWIDTH,
            linecolor='black',
            )
# plt.show()
plt.savefig(MODEL+'/'+INPUT_FILE.split('.')[0]+'_mutype.svg', dpi=1000)