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
from tqdm import tqdm
sys.path.insert(1, '../ML/')
import fetchData

import numpy as np
from matplotlib_venn import venn3, venn3_circles

class HMM:
	def __init__(self, pfampos) -> None:
		self.pfampos = pfampos
		self.mut_types = {}
		self.ptm_types = 0
		self.bitscore = None

class Mutation:
	def __init__(self, acc, mutation, pfampos, mut_type) -> None:
		self.mutation = mutation
		self.position = int(mutation[1:-1])
		self.acc = acc
		self.pfampos = pfampos
		self.mut_type = mut_type
		self.ptm_types = []
		self.bitscore = None
		self.orthscore = None

mydb = fetchData.connection(db_name='kinase_project2')
mydb.autocommit = True
mycursor = mydb.cursor()

dic_hmm = {}
dic_mutations = {}
'''Load mutation information for a pfam_pos'''
mycursor.execute("SELECT acc, mutation, pfampos, mut_type FROM mutations\
		 		where pfampos!=%s and pfampos!=%s", ('-', 'neutral'))
mutations = mycursor.fetchall()
inc = 0
dec = 0
gain = 0
loss = 0
res = 0
for mut_row in mutations:
	acc, mutation, pfampos, mut_type = mut_row
	name = acc + '_' + mutation
	pfampos = int(pfampos)
	if mut_type == 'increase':
		inc += 1
	elif mut_type == 'decrease':
		dec += 1
	elif mut_type == 'activating':
		gain += 1
	elif mut_type == 'loss':
		loss += 1
	elif mut_type == 'resistance':
		res +=1

data = []
data.append([gain, 'constitutive\nactivation', 'activating'])
data.append([inc, 'increase', 'activating'])
data.append([loss, 'loss', 'deactivating'])
data.append([dec, 'decrease', 'deactivating'])
data.append([res, 'resistance', 'resistance'])

df = pd.DataFrame(data, columns=['# mutations', 'Mutation type', 'mut_type2'])
ax = sns.barplot(data=df, x="Mutation type", y="# mutations", palette=['green', 'lightgreen', 'red', 'lightcoral', 'cornflowerblue'])
for i in ax.containers:
    ax.bar_label(i,)
plt.grid(linewidth=0.5, color='gray', linestyle='--')
# ax.legend(loc="upper center")
plt.savefig('mutation_count.png', dpi=1000)
plt.savefig('mutation_count.eps', dpi=1000)
plt.savefig('mutation_count.svg', dpi=1000)
# plt.show()