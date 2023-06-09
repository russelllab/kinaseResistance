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
dic_mutations = {'upstream AL': {}, 'AL': {}, 'downstream AL': {}}
'''Load mutation information for a pfam_pos'''
mycursor.execute("SELECT acc, mutation, pfampos, mut_type FROM mutations\
		 		where pfampos!=%s and pfampos!=%s", ('-', 'neutral'))
mutations = mycursor.fetchall()
for mut_row in mutations:
	acc, mutation, pfampos, mut_type = mut_row
	if mut_type == 'activating': mut_type = 'constitutive-activation'
	if pfampos == '-': continue
	name = acc + '_' + mutation
	pfampos = int(pfampos)
	if pfampos < 442:
		if mut_type not in dic_mutations['upstream AL']: dic_mutations['upstream AL'][mut_type] = 0
		dic_mutations['upstream AL'][mut_type] += 1
	elif 442 <= pfampos <= 534:
		if mut_type not in dic_mutations['AL']: dic_mutations['AL'][mut_type] = 0
		dic_mutations['AL'][mut_type] += 1
	else:
		if mut_type not in dic_mutations['downstream AL']: dic_mutations['downstream AL'][mut_type] = 0
		dic_mutations['downstream AL'][mut_type] += 1


pfam_positions = [pfam_pos for pfam_pos in list(set(dic_hmm.keys()))]
pfam_positions.sort()

fig, ax = plt.subplots()

mut_types = ['constitutive-activation', 'increase', 'resistance', 'decrease', 'loss']
mut_types_colors = ['green', 'lightgreen', 'blue', 'lightcoral', 'red']
bottom = np.zeros(3)
for mut_type, mut_types_color in zip(mut_types, mut_types_colors):
	row = []
	for region in ['upstream AL', 'AL', 'downstream AL']:
		num = dic_mutations[region][mut_type] if mut_type in dic_mutations[region] else 0
		row.append(num)
	row = np.array(row)
	print (mut_type, row)
	width = 0.5
	p = ax.bar(['N-lobe', 'AL', 'C-lobe'], row, width,\
	    		label=mut_type, bottom=bottom, color=mut_types_color)
	bottom += row

ax.set_title("Distribution of mutations")
ax.legend(loc="upper center")
# plt.show()
plt.grid(linewidth=0.5, color='gray', linestyle='--')
plt.savefig('distribution_of_mutations_KD.png', dpi=1000)
plt.savefig('distribution_of_mutations_KD.eps', dpi=1000)
plt.savefig('distribution_of_mutations_KD.svg', dpi=1000)