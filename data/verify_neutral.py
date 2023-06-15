#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
A script to create a circular bar plot.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
sys.path.append('../ML/')
import fetchData
from tqdm import tqdm
import seaborn as sns

AA = 'ACDEFGHIKLMNPQRSTVWY'
mydb = fetchData.connection(db_name='kinase_project2')
mydb.autocommit = True
mycursor = mydb.cursor()

dic_ptms = {}
mycursor.execute("SELECT uniprotpos, gene, ptmtype FROM ptms")
hits = mycursor.fetchall()
for hit in hits:
    uniprotpos, gene, ptmtype = hit
    uniprotpos = int(uniprotpos)
    if gene not in dic_ptms:
        dic_ptms[gene] = {}
    if uniprotpos not in dic_ptms[gene]:
        dic_ptms[gene][uniprotpos] = []
    dic_ptms[gene][uniprotpos].append(ptmtype)

mycursor.execute("SELECT mutation, gene, acc, mut_type FROM mutations\
                 where pfampos!='%s'" % ('-', ))
hits = mycursor.fetchall()
dic_var = {}
for hit in hits:
    mutation, gene, acc, mut_type = hit
    if mut_type not in dic_var:
        dic_var[mut_type] = []
    dic_var[mut_type].append(acc+'/'+mutation)
    if mut_type == 'neutral':
        wtPos = int(mutation[1:-1])
        if gene in dic_ptms:
            if wtPos in dic_ptms[gene]:
                print (f'Neutral mutation {mutation} is a {dic_ptms[gene][wtPos]}-site')

for var in dic_var['neutral']:
    for mut_type in dic_var:
        if mut_type == 'neutral': continue
        for var2 in dic_var[mut_type]:
            if var == var2:
                print (f'Neutral mutation {var} is also a {mut_type} mutation')

df_dic = {}
for mut_type in dic_var:
    if mut_type not in df_dic:
        df_dic[mut_type] = []
    for var in tqdm(dic_var[mut_type]):
        acc = var.split('/')[0]
        mutation = var.split('/')[1]        
        wtPos = int(mutation[1:-1])
        wtAA = mutation[0]
        mutAA = mutation[-1]

        mycursor.execute("SELECT * FROM spec_para where acc=%s and position=%s",\
                         (acc, wtPos))
        hit = mycursor.fetchone()
        for aa, value in zip(AA, hit[3:-1]):
            if aa == wtAA:
                wtValue = float(value)
            elif aa == mutAA:
                mutValue = float(value)
        diff = mutValue - wtValue
        df_dic[mut_type].append(diff)

data = []
for mut_type in df_dic:
    for diff in df_dic[mut_type]:
        data.append([mut_type, diff])

df = pd.DataFrame(data, columns=['mut_type', 'diff'])
print (df)

sns.violinplot(x='mut_type', y='diff', data=df)
plt.grid(linewidth=0.5)
plt.show()