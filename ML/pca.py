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

df = pd.read_csv('trainData.tsv.gz', sep = '\t')
# exclude columns
# df = df.loc[:, ~df.columns.isin(['allHomologs','exclParalogs','specParalogs','orthologs', 'bpso','bpsh'])]
df = df.loc[:, ~df.columns.isin([
                            'allHomologs',
                            'exclParalogs',
                            'specParalogs',
                            # 'orthologs',
                            'bpso',
                            'bpsh'
                            ])]
# exclude columns to make the data matrix
original_df = df.copy()
df = df.loc[:, ~df.columns.isin(['Acc', 'Mutation'])]
print (df.to_numpy())
# sys.exit()
data = []
mut_types_colors = []
# for line in gzip.open('trainData.tsv.gz', 'rt'):
for row in df.to_numpy():
    mut_type = row[-1]
    if mut_type == 'A':
        # continue
        mut_types_colors.append('green')
    elif mut_type == 'D':
        mut_types_colors.append('red')
    else:
        # continue
        mut_types_colors.append('violet')
    data.append(row[:-1])

# data = np.array(data)
# scaler = MinMaxScaler()
# scaler.fit(data)
# data = scaler.transform(data)
# print (data[:,0])
# print (trainMat)
# sys.exit()

pca = decomposition.PCA(n_components=2)
pca.fit(data)
X = pca.transform(data)
print (pca.explained_variance_ratio_, pca.n_components_)

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

pca_df = pd.DataFrame(X, columns = ['PCA1', 'PCA2'])
pca_df = pd.concat([pca_df, original_df], axis=1)
print (pca_df)
# print (original_df)
plt.show()

fig = px.scatter(
                pca_df, x="PCA1", y="PCA2", color="MUT_TYPE",
                 hover_data=['Acc', 'Mutation', 'p', 'p_pfam'],
                 color_discrete_map={
                                    "A": "green",
                                    "D": "red",
                                    "R": "blue"
                                    }
                 )
fig.show()
# plt.savefig('pca_plot.png')