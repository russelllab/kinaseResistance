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

mycursor.execute("SELECT pfampos, pfamaa from hmm")
pfam = mycursor.fetchall()
dic_pfam = {}
for pfam_row in pfam:
    pfampos, pfamaa = pfam_row
    if pfampos not in dic_pfam: dic_pfam[pfampos] = pfamaa

# extracts p-sites
mycursor.execute("SELECT acc, uniprotpos, ptmtype FROM ptms")
ptms = mycursor.fetchall()
dic_ptms = {}
for ptm_row in ptms:
    acc, uniprotpos, ptmtype = ptm_row
    uniprotpos = int(uniprotpos)
    if ptmtype != 'p': continue
    if acc not in dic_ptms: dic_ptms[acc] = []
    dic_ptms[acc].append(uniprotpos)

# extracts mutations
mycursor.execute("SELECT acc, wtpos, mutation, pfampos, mut_type FROM mutations\
                    where mut_type!=%s", ('neutral',))
mutations = mycursor.fetchall()
dic_hmm = {}
for mutation in mutations:
    acc, wtpos, mutation, pfampos, mut_type = mutation
    if mut_type =='activating': mut_type = 'constitutive-activation'
    if pfampos == '-': continue
    mutAA = mutation[-1]
    print (mutAA)
    if acc in dic_ptms:
        if int(wtpos) in dic_ptms[acc]:
            mut_type += ' at p-sites'
            if mutAA in ['D', 'E']:
                mut_type += ' to DE'
    if pfampos not in dic_hmm: dic_hmm[pfampos] = {}
    if mut_type not in dic_hmm[pfampos]: dic_hmm[pfampos][mut_type] = 0
    dic_hmm[pfampos][mut_type] += 1

data = []
for hmmpos in dic_hmm:
    total = 0
    for mut_type in dic_hmm[hmmpos]:
        total += dic_hmm[hmmpos][mut_type]
    if total >= 10:
        # print(hmmpos, total, dic_hmm[hmmpos])
        total_log2 = np.log2(total)
        for mut_type in dic_hmm[hmmpos]:
            mut_type_count = dic_hmm[hmmpos][mut_type]
            mut_type_fraction = float(mut_type_count) / total
            mut_type_count_log2 = total_log2 * mut_type_fraction
            dic_hmm[hmmpos][mut_type] = mut_type_count_log2
        pfamaa = dic_pfam[hmmpos]
        data.append([hmmpos, pfamaa, total, dic_hmm[hmmpos]])

df = pd.DataFrame(data, columns=['pfampos', 'pfamaa', 'total', 'mut_types'])
df = df.sort_values(by=['total'], ascending=False)
print (df)
print (df['pfampos'].to_numpy())
pfam_positions = df['pfampos'].to_numpy()
pfam_position_labels = []
for pfam_position in pfam_positions:
    pfam_position_labels.append(str(pfam_position) + ':' + dic_pfam[pfam_position])
mut_types = ['constitutive-activation', 'increase', 'resistance', 'decrease', 'loss']
new_mut_types = []
for mut_type in mut_types:
    new_mut_types.append(mut_type)
    new_mut_types.append(mut_type+' at p-sites')
    new_mut_types.append(mut_type+' at p-sites to DE')
mut_types = new_mut_types
mut_types_colors = ['green', 'lightgreen', 'blue', 'lightcoral', 'red']
new_mut_types_colors = []
for mut_type_color in mut_types_colors:
    new_mut_types_colors.append(mut_type_color)
    new_mut_types_colors.append(mut_type_color)
    new_mut_types_colors.append(mut_type_color)
mut_types_colors = new_mut_types_colors
left = np.zeros(len(pfam_positions))
width = 0.5
fig, ax = plt.subplots()
hatches = []
for mut_type in mut_types:
    print (mut_type)
    # hatch='.....' if ' at p-sites' in mut_type else None
    if ' at p-sites to DE' in mut_type:
        hatches.append('///.....')
    elif ' at p-sites' in mut_type:
        hatches.append('.....')
        # hatches.append(None)
    else:
        hatches.append(None)
print (hatches)
for mut_type, mut_type_color, hatch in zip(mut_types, mut_types_colors, hatches):
    row = []
    for pfam_position in pfam_positions:
        if mut_type not in dic_hmm[pfam_position]:
            row.append(0)
        else:
            row.append(dic_hmm[pfam_position][mut_type])
    row = np.array(row)
    p = ax.barh(pfam_position_labels, row, width,\
                label=mut_type, left=left, color=mut_type_color,\
                edgecolor = 'black', linewidth = 0.5, \
                hatch = hatch
                # hatch='.....' if ' at p-sites' in mut_type else None
                )
    left += row

ax.set_title("Top 15 most mutated alignment positions")
ax.set_xlim(0,8)
# ax.legend(loc="upper right")
plt.grid(linewidth=0.5, color='gray', linestyle='--')
# plt.show()
plt.savefig('most_mutated_pfam_positions.png', dpi=1000)
plt.savefig('most_mutated_pfam_positions.eps', dpi=1000)
plt.savefig('most_mutated_pfam_positions.svg', dpi=1000)