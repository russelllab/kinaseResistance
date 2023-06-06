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
sys.path.insert(1, '../ML/')
import fetchData

class HMM:
	def __init__(self, pfampos) -> None:
		self.pfampos = pfampos
		self.mut_types = {}
		self.ptm_types = 0
		self.bitscore = None

mydb = fetchData.connection(db_name='kinase_project2')
mydb.autocommit = True
mycursor = mydb.cursor()

'''Load ptm information for a pfam_pos'''
dic_ptm_pfam = {}
mycursor.execute("SELECT pfampos, ptmtype FROM ptms\
		 		where pfampos!='-'")
ptms = mycursor.fetchall()
dic_ptm_pfam = {}
dic_hmm = {}
for ptm_row in ptms:
	pfampos = int(ptm_row[0])
	ptm_type = ptm_row[1]
	if pfampos not in dic_hmm: dic_hmm[pfampos] = HMM(pfampos)
	dic_hmm[pfampos].ptm_types += 1

mycursor.execute("SELECT mutation, pfampos, mut_type FROM mutations\
		 		where pfampos!=%s and pfampos!=%s", ('-', 'neutral'))
mutations = mycursor.fetchall()
# dic_mut_types = {'increase': 'A'
for mut_row in mutations:
	mutation, pfampos, mut_type = mut_row
	pfampos = int(pfampos)
	if pfampos not in dic_hmm: dic_hmm[pfampos] = HMM(pfampos)
	if mut_type not in dic_hmm[pfampos].mut_types:
		dic_hmm[pfampos].mut_types[mut_type] = []
	dic_hmm[pfampos].mut_types[mut_type].append(mutation)


'''Load HMM information for a pfam_pos'''
mycursor.execute("SELECT * FROM hmm")
hmm = mycursor.fetchall()
for hmm_row in hmm:
	pfampos, pfamaa = hmm_row[0], hmm_row[1]
	if pfampos == '-': continue
	pfampos = int(pfampos)
	if pfampos not in dic_hmm: continue
	num = 4
	for count, aa in enumerate(['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
				'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']):
		if aa == pfamaa.upper():
			bitscore = float(hmm_row[num+count])
			dic_hmm[pfampos].bitscore = bitscore
			break

# Normalize data
for pfampos in dic_hmm:
	# numpy.log2
	for mut_type in dic_hmm[pfampos].mut_types:
		# print (dic_hmm[pfampos].mut_types[mut_type])
		if len(dic_hmm[pfampos].mut_types[mut_type]) <2:
			dic_hmm[pfampos].mut_types[mut_type] = 1
		else:
			dic_hmm[pfampos].mut_types[mut_type] = np.log2(len(dic_hmm[pfampos].mut_types[mut_type]))
	'''
	total = 0
	for mut_type in dic_hmm[pfampos].mut_types:
		total += len(dic_hmm[pfampos].mut_types[mut_type])
	for mut_type in dic_hmm[pfampos].mut_types:
		# print (dic_hmm[pfampos].mut_types[mut_type])
		dic_hmm[pfampos].mut_types[mut_type] = len(dic_hmm[pfampos].mut_types[mut_type]) / float(total)
		dic_hmm[pfampos].mut_types[mut_type] = round(dic_hmm[pfampos].mut_types[mut_type], 2)
	'''

pfam_positions = [pfam_pos for pfam_pos in list(set(dic_hmm.keys()))]

col_colors = []
for pfampos in pfam_positions:
	# print (pfampos)
	if dic_hmm[pfampos].bitscore <= 0.5:
		col_colors.append('red')
		print (pfampos, dic_hmm[pfampos].bitscore)
	elif dic_hmm[pfampos].ptm_types >= 10:
		col_colors.append('grey')
	else:
		col_colors.append('white')

mut_types = ['activating', 'increase', 'resistance', 'loss', 'decrease']
data = []
for mut_type in mut_types:
	row = []
	for pfam_pos in pfam_positions:
		if mut_type not in dic_hmm[pfam_pos].mut_types:
			value = 0
		else:
			value = dic_hmm[pfam_pos].mut_types[mut_type]
		row.append(value)
	data.append(row)

# print (data)
# print ([pfam_pos for pfam_pos in pfam_positions])

df = pd.DataFrame(
				data,
				index = mut_types,
				columns = pfam_positions
				)
# df = df.loc[:, (df != 0).any(axis=0)]

# print (df)

palette = sns.light_palette("seagreen", as_cmap=True)
g = sns.clustermap(df, cmap=palette, row_cluster=False, col_cluster=False,\
		   			col_colors=col_colors)
#plt.savefig('clusterMap.jpeg', format="jpeg", dpi=500)
plt.show()
sys.exit()