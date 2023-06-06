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

os.system('hmmsearch -o out.txt ../pfam/Pkinase.hmm humanKinases.fasta')
ACC = 'Q9NYV4'
flag = 0
dic_seq2pfam = {}
for line in open('out.txt', 'r'):
    if line.startswith('>>'):
        acc = line.split()[1].split('|')[1]
        if acc == ACC:
            print (acc)
            flag = 1
        else:
            flag = 0
    elif flag == 1:
        if line.strip() == '': continue
        if line.split()[0] == 'Pkinase':
            hmm_start = int(line.split()[1])
            hmm_end = int(line.split()[-1])
            hmm_seq = line.split()[2]
        elif acc in line.split()[0]:
            seq_start = int(line.split()[1])
            seq_end = int(line.split()[-1])
            seq_seq = line.split()[2]
            for hmm_res, seq_res in zip(hmm_seq, seq_seq):
                if hmm_res in ['.', '-']:
                    seq_start += 1
                    continue
                if seq_res in ['.', '-']:
                    hmm_start += 1
                    continue
                # print (hmm_res, seq_res, hmm_start, seq_start)
                dic_seq2pfam[seq_start] = hmm_start
                hmm_start += 1
                seq_start += 1

mycursor.execute("SELECT uniprotaa, uniprotpos, pfampos, acc from positions\
    where acc='%s'" % (ACC))
positions = mycursor.fetchall()
for position in positions:
    uniprotaa, uniprotpos, pfampos, acc = position
    uniprotpos = int(uniprotpos)
    if uniprotpos in dic_seq2pfam:
        print (uniprotaa, uniprotpos, pfampos, dic_seq2pfam[uniprotpos])