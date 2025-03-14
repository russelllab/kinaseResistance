#!/usr/bin/env python3
# -*- coding: utf-8 -*-

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
from sklearn.manifold import TSNE

AA = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

df = pd.read_csv('trainDataFromHitsSplitTrimmedAln.tsv.gz', sep = '\t')

# Enable this if you want to plot only train data
# df = df[df.Dataset == 'train']
# print (df.to_numpy().tolist())
# sys.exit()

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
                    'Acc',
                    'Mutation',
                    'Gene',
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
                    #   'Acetylmimic',
                    #   'hmmScoreWT',
                    #   'hmmScoreMUT',
                    #   'hmmScoreDiff'
                    ]

columns_to_exclude += ['ncontacts', 'nresidues', 'mech_intra']
columns_to_exclude += ['phi_psi', 'sec', 'burr', 'acc']
columns_to_exclude += ['IUPRED']

for aa in AA:
    # if aa not in ['S', 'T', 'Y']:
        columns_to_exclude.append(aa+'_WT')
    # if aa not in ['D', 'E']:
        columns_to_exclude.append(aa+'_MUT')


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
    # if i in [-1, 0, 1]: continue
    for col in adr_cols:
        columns_to_exclude.append(col+'_'+str(i))

############

adr_cols = ['A_pfam', 'D_pfam', 'R_pfam']
# adr_cols = ['D_pfam', 'R_pfam']
for i in range(-5, 6):
    # if i in [-1, 0, 1]: continue
    for col in adr_cols:
        columns_to_exclude.append(col.split('_')[0]+'_'+str(i)+'_'+col.split('_')[1])

adr_cols = ['A', 'D', 'R']
AA = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 
    'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
for mut_type in adr_cols:
    for aa in AA:
        col = mut_type+'_'+aa
        for i in range(-5, 6):
            if i in [0]: continue
            columns_to_exclude.append(col+'_'+str(i)+'_pfam')


df = df.loc[:, ~df.columns.isin(columns_to_exclude)]
# print (df)

# scaler = MinMaxScaler()
# columns_to_scale = ['p_pfam', 'ac_pfam', 'me_pfam', 'gl_pfam', 'm1_pfam', 'm2_pfam', 'm3_pfam', 'sm_pfam', 'ub_pfam']
# columns_to_scale += ['hmmScoreDiff', 'hmmScoreWT', 'hmmScoreMUT']
# df[columns_to_scale] = scaler.fit_transform(df[columns_to_scale])

print (df)
# sys.exit()
data = []
mut_types_colors = []
# for line in gzip.open('trainData.tsv.gz', 'rt'):
for row in df.to_numpy():
    mut_type = row[-1]
    if mut_type in ['activating', 'increase']:
        # continue
        mut_types_colors.append('green')
    elif mut_type in ['loss', 'decrease']:
        mut_types_colors.append('red')
    elif mut_type == 'resistance':
        mut_types_colors.append('blue')
    elif mut_type == 'neutral':
        mut_types_colors.append('cyan')
    elif mut_type == 'activatingresistance':
        mut_types_colors.append('violet')
    else:
        # continue
        mut_types_colors.append('yellow')
    data.append(row[:-1])

data = np.array(data)
scaler = MinMaxScaler()
scaler.fit(data)
data = scaler.transform(data)
# print (data[:,0])
# print (trainMat)
# sys.exit()

X = TSNE(n_components=2, learning_rate='auto', init='random', perplexity=100).fit_transform(data)

# pca = decomposition.PCA(n_components=2)
# pca.fit(data)
# X = pca.transform(data)
# print (pca.explained_variance_ratio_, pca.n_components_)

fig = plt.figure(1, figsize=(4, 3))
plt.clf()

ax = fig.add_subplot(111)
# ax = fig.add_subplot(111, projection="3d", elev=48, azim=134)
ax.set_position([0.1, 0.1, 0.8, 0.8])


plt.cla()

# ax.scatter(X[:, 0], X[:, 1], X[:, 2], c=mut_type_colors, cmap=plt.cm.nipy_spectral, edgecolor="k")
ax.scatter(X[:, 0], X[:, 1], c=mut_types_colors, cmap=plt.cm.nipy_spectral, edgecolor="k")

kmeans = KMeans(n_clusters=2, random_state=0).fit(X)
print (kmeans.labels_)
# kmeans.predict([[0, 0], [12, 3]])
print (kmeans.cluster_centers_)

col1 = 'PCA1'
col1 = 'tsne_1'
col2 = 'PCA2'
col2 = 'tsne_2'
pca_df = pd.DataFrame(X, columns = [col1, col2])
pca_df = pd.concat([pca_df, original_df], axis=1)
# print (pca_df)
# print (original_df)
plt.show()
# plt.savefig('pca_plot.png')

pca_df.to_csv('trainDataPostPCA.tsv.gz', sep = '\t')

fig = px.scatter(
                pca_df, x="tsne_1", y="tsne_2", color="MUT_TYPE",
                # symbol = 'Dataset',
                # size = pca_df['Dataset'].to_numpy(),
                size = 'Dataset',
                hover_data=[
                        'Gene',
                        'Mutation',
                        # 'Phosphomimic',
                        'p_0',
                        # 'p_0_pfam',
                        # 'ac',
                        # 'ac_pfam',
                        # 'me',
                        # 'me_pfam',
                        # 'hmmScoreWT',
                        # 'hmmScoreMUT',
                        # 'hmmSS',
                        'hmmPos'
                        ],
                color_discrete_map={
                                "activating": "green",
                                "increase": "lightgreen",
                                "loss": "red",
                                "decrease": "lightcoral",
                                "resistance": "blue",
                                "neutral": "cyan",
                                "activatingresistance": "violet",
                                "TBD": "yellow",
                                "A": "darkgreen",
                                "Inconclusive": "grey"
                                }
                 )
fig.show()