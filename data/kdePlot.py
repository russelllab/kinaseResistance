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
import plotly.graph_objects as go
sys.path.insert(1, '../ML/')
import fetchData

import numpy as np
from matplotlib_venn import venn3, venn3_circles

mydb = fetchData.connection(db_name='kinase_project2')
mydb.autocommit = True
mycursor = mydb.cursor()

mycursor.execute("SELECT mut_type, pfampos, acc from mutations\
                 where pfampos!=%s", ('-',))
mutations = mycursor.fetchall()
dic_hmm = {}
for mutation in mutations:
    mut_type, pfampos, acc = mutation
    if mut_type == 'neutral': continue
    if mut_type in ['activating', 'increase']:
        mut_type_group = 'activating'
    elif mut_type == 'loss' or mut_type == 'decrease':
        mut_type_group = 'deactivating'
    elif mut_type == 'resistance':
        mut_type_group = 'resistance'
    pfampos = int(pfampos)
    if pfampos not in dic_hmm: dic_hmm[pfampos] = {}
    if mut_type_group not in dic_hmm[pfampos]: dic_hmm[pfampos][mut_type_group] = 0
    dic_hmm[pfampos][mut_type_group] += 1

print (dic_hmm)

data = []
for pfampos in dic_hmm:
    for mut_type_group in dic_hmm[pfampos]:
        data.append([pfampos, mut_type_group, dic_hmm[pfampos][mut_type_group]])
df = pd.DataFrame(data, columns=['pfampos', 'mut_type_group', 'count'])

sns.histplot(
data=df, x="pfampos", hue="mut_type_group",
fill=True, common_norm=False,
palette={'activating':'green', 'deactivating':'red', 'resistance':'cornflowerblue'},
element='step', bins=200, multiple='stack',
# clip=(30,764), bw_adjust=0.05,

alpha=.5, linewidth=0,
)
# plt.savefig('kdeplot_'+mut_type_group+'.png', dpi=1000)
plt.show()
# ax.clear()
# plt.clf()
# plt.close()
