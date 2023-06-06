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

mycursor.execute("SELECT uniprotpos, pfampos, acc from positions")
positions = mycursor.fetchall()
dic_pos = {}
for position in positions:
    uniprotpos, pfampos, acc = position
    if pfampos == '-': continue
    uniprotpos = int(uniprotpos)
    pfampos = int(pfampos)
    if acc not in dic_pos: dic_pos[acc] = {'start': None, 'end': None}
    if dic_pos[acc]['start'] == None:
        dic_pos[acc]['start'] = uniprotpos
    else:
        if uniprotpos < dic_pos[acc]['start']: dic_pos[acc]['start'] = uniprotpos
    if dic_pos[acc]['end'] == None:
        dic_pos[acc]['end'] = uniprotpos
    else:
        if uniprotpos > dic_pos[acc]['end']: dic_pos[acc]['end'] = uniprotpos

# print (dic_pos)


mycursor.execute("SELECT acc, mutation, wtpos, pfampos, mut_type FROM mutations\
                    where mut_type!=%s", ('neutral',))
mutations = mycursor.fetchall()
dic_kd = {}
for mutation in mutations:
    acc, mutation, wtpos, pfampos, mut_type = mutation
    wtpos = int(wtpos)
    if mut_type not in dic_kd: dic_kd[mut_type] = {'kd':0, 'ukd':0, 'dkd':0}
    if pfampos != '-':
        print (pfampos)
        dic_kd[mut_type]['kd'] += 1
    else:
        if acc not in dic_pos: continue
        if wtpos < dic_pos[acc]['start']:
            dic_kd[mut_type]['ukd'] += 1
        elif wtpos > dic_pos[acc]['end']:
            dic_kd[mut_type]['dkd'] += 1

print (dic_kd)
# sys.exit()
level1 = ['COSMIC', 'UniProt', 'Text Mining']
level2 = ['Constitutive-activation', 'Increase', 'Loss', 'Decrease', 'Resistance']
level3 = ['Upstream KD', 'KD', 'Downstream KD']
levels = level1 + level2 + level3
level_cosmic_resistance = 98
level_uniprot_activating = 129
level_uniprot_increase = 203
level_uniprot_loss = 567
level_uniprot_decrease = 468
level_text_mining_activating = 38
# print (len(level_cosmic_resistance))

fig = go.Figure(data=[go.Sankey(
    node = dict(
      pad = 15,
      thickness = 20,
      line = dict(color = "black", width = 0.5),
    #   label = ["A1", "A2", "B1", "B2", "C1", "C2"],
      label = levels,
      color = ['gold', 'cyan', 'grey',\
                'green','lightgreen', 'red', 'lightcoral','cornflowerblue',\
                'grey', 'grey', 'grey']
    ),
    link = dict(
      source = [0, 1, 1, 1, 1, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7], # indices correspond to labels, eg A1, A2, A1, B1, ...
      target = [7, 3, 4, 5, 6, 3, 8, 9, 10, 8, 9, 10, 8, 9, 10, 8, 9, 10, 8, 9],
      value = [98, 129, 203, 567, 468, 38, 45, 99, 18, 43, 93,48,61, 404, 34, 120, 232, 75, 10, 86],
      color = ['gold', 'cyan', 'cyan', 'cyan', 'cyan', 'silver', 'green', 'green', 'green'\
               ,'lightgreen', 'lightgreen', 'lightgreen', 'red', 'red', 'red'\
                ,'lightcoral','lightcoral','lightcoral'\
                ,'cornflowerblue','cornflowerblue','cornflowerblue'\
                    ]
  ))])

fig.update_layout(title_text="Basic Sankey Diagram", font_size=10)
fig.show()