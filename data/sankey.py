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

mycursor.execute("SELECT uniprotpos, alnpos, pfampos, acc from positions")
positions = mycursor.fetchall()
dic_pos = {}
for position in positions:
    uniprotpos, alnpos, pfampos, acc = position
    if pfampos == '-': continue
    # alnpos = int(alnpos)
    alnpos = int(pfampos)
    # if alnpos not in range(30, 833): continue
    uniprotpos = int(uniprotpos)
    pfampos = int(pfampos)
    if acc not in dic_pos: dic_pos[acc] = {'ntkdstart': None,'kdstart': None, 'kdend': None, 'ctkdend': None}
    if alnpos in range(30, 833):
        if dic_pos[acc]['kdstart'] == None:
            dic_pos[acc]['kdstart'] = uniprotpos
        else:
            if uniprotpos < dic_pos[acc]['kdstart']: dic_pos[acc]['kdstart'] = uniprotpos
        if dic_pos[acc]['kdend'] == None:
            dic_pos[acc]['kdend'] = uniprotpos
        else:
            if uniprotpos > dic_pos[acc]['kdend']: dic_pos[acc]['kdend'] = uniprotpos
    else:
        if alnpos < 30:
            if dic_pos[acc]['ntkdstart'] == None:
                dic_pos[acc]['ntkdstart'] = uniprotpos
            else:
                if uniprotpos < dic_pos[acc]['ntkdstart']: dic_pos[acc]['ntkdstart'] = uniprotpos
        else:
            if dic_pos[acc]['ctkdend'] == None:
                dic_pos[acc]['ctkdend'] = uniprotpos
            else:
                if uniprotpos > dic_pos[acc]['ctkdend']: dic_pos[acc]['ctkdend'] = uniprotpos

# print (dic_pos)
# sys.exit()
for acc in dic_pos:
    if dic_pos[acc]['ntkdstart'] == None: dic_pos[acc]['ntkdstart'] = dic_pos[acc]['kdstart']
    if dic_pos[acc]['ctkdend'] == None: dic_pos[acc]['ctkdend'] = dic_pos[acc]['kdend']

mycursor.execute("SELECT acc, mutation, wtpos, pfampos, mut_type, source FROM mutations\
                    where mut_type!=%s", ('neutral',))
mutations = mycursor.fetchall()
dic_kd = {}
dic_source = {}
for mutation in mutations:
    acc, mutation, wtpos, pfampos, mut_type, source = mutation
    wtpos = int(wtpos)
    if mut_type not in dic_kd: dic_kd[mut_type] = {'kd':0, 'nkd':0, 'ckd':0, 'ntkd':0, 'ctkd':0}
    # if pfampos != '-':
    #     print (pfampos)
    #     dic_kd[mut_type]['kd'] += 1
    # else:
    if acc not in dic_pos: continue

    if source not in dic_source: dic_source[source] = {}
    if mut_type not in dic_source[source]: dic_source[source][mut_type] = 0
    dic_source[source][mut_type] += 1

    print (acc, wtpos, dic_pos[acc]['kdstart'])
    if wtpos < dic_pos[acc]['ntkdstart']:
        dic_kd[mut_type]['nkd'] += 1
    elif wtpos >= dic_pos[acc]['ntkdstart'] and wtpos < dic_pos[acc]['kdstart']:
        dic_kd[mut_type]['ntkd'] += 1
    elif wtpos >= dic_pos[acc]['kdstart'] and wtpos <= dic_pos[acc]['kdend']:
        dic_kd[mut_type]['kd'] += 1
    elif wtpos > dic_pos[acc]['kdend'] and wtpos <= dic_pos[acc]['ctkdend']:
        dic_kd[mut_type]['ctkd'] += 1
    else:
        dic_kd[mut_type]['ckd'] += 1

print (dic_kd)
print (dic_source)
# sys.exit()
level1 = ['COSMIC', 'UniProt', 'PubMed', 'UniProt+PubMed']
level2 = ['activating', 'increase', 'loss', 'decrease', 'resistance']
level3 = ['nkd', 'ntkd', 'kd', 'ctkd', 'ckd']
levels = level1 + level2 + level3
dic_colors = {'COSMIC': 'gold', 'UniProt': 'lightblue', 'PubMed': 'silver',
              'UniProt+PubMed': 'grey', 'activating': 'green', 'increase': 'lightgreen',
              'loss': 'red', 'decrease': 'lightcoral', 'resistance': 'cornflowerblue'}
level_cosmic_resistance = 98
level_uniprot_activating = 129
level_uniprot_increase = 203
level_uniprot_loss = 567
level_uniprot_decrease = 468
level_text_mining_activating = 38
# print (len(level_cosmic_resistance))
source = []
target = []
value = []
color = []
for l1 in level1:
    for l2 in level2:
        if l2 not in dic_source[l1]: continue
        source.append(levels.index(l1))
        target.append(levels.index(l2))
        value.append(dic_source[l1][l2])
        color.append(dic_colors[l1])

for l2 in level2:
    for l3 in level3:
        source.append(levels.index(l2))
        target.append(levels.index(l3))
        value.append(dic_kd[l2][l3])
        color.append(dic_colors[l2])

print (source)
print (target)
print (value)
print (color)
# sys.exit()


fig = go.Figure(data=[go.Sankey(
    node = dict(
      pad = 5,
      thickness = 10,
      line = dict(color = "black", width = 0.5),
    #   label = ["A1", "A2", "B1", "B2", "C1", "C2"],
      label = levels,
      x = [0.05, 0.05, 0.05, 0.05,
           0.5, 0.5, 0.5, 0.5, 0.5,
           1.0, 1.0, 1.0, 1.0, 1.0],
      y = [1.1, 0.45, 1.0, 0.94,
           0.8, 0.95, 0.22, 0.57, 1.075,
           0.1, 0.25, 0.64, 1.01, 1.1],
      color = ['gold', 'lightblue', 'grey', 'grey',
                'green','lightgreen', 'red', 'lightcoral','cornflowerblue',
                'grey', 'grey', 'grey', 'grey', 'grey']
    ),
    link = dict(
        source = source,
        target = target,
        value =  value,
        color = color
    )
)])
'''
    link = dict(
      source = [0,    1,   1,   1,   1,  2,  3,  4,  4,  4,  5,  5,  5,  6,   6,  6,   7,   7,  7,  8,  8], # indices correspond to labels, eg A1, A2, A1, B1, ...
      target = [8,    4,   5,   6,   7,  4,  4,  9, 10, 11,  9, 10, 11,  9,   10, 11,  9,  10, 11,  9, 10],
      value =  [98, 129, 184, 499, 427, 38,  9, 45, 99, 18, 43, 93, 48, 61, 404, 34, 120, 232, 75, 10, 86],
      color = ['gold', 'lightblue', 'lightblue', 'lightblue', 'lightblue', 'silver', 'silver', 'green', 'green', 'green'\
               ,'lightgreen', 'lightgreen', 'lightgreen', 'red', 'red', 'red'\
                ,'lightcoral','lightcoral','lightcoral'\
                ,'cornflowerblue','cornflowerblue','cornflowerblue'\
                    ]
    )
    
  )])
'''

fig.update_layout(title_text="Basic Sankey Diagram", font_size=10)
# fig.show()
fig.write_image("sankey_plot.svg")