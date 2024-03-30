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
		self.neighbour_pac_sites = 0
		self.bitscore = None
		self.orthscore = None

mydb = fetchData.connection(db_name='kinase_project2')
mydb.autocommit = True
mycursor = mydb.cursor()

dic_hmm = {}
dic_mutations = {}
'''Load mutation information for a pfam_pos'''
mycursor.execute("SELECT acc, mutation, pfampos, mut_type FROM mutations\
		 		where pfampos!=%s", ('-', ))
mutations = mycursor.fetchall()
for mut_row in mutations:
	acc, mutation, pfampos, mut_type = mut_row
	if mut_type =='activating': mut_type = 'constitutive-activation'
	name = acc + '_' + mutation
	pfampos = int(pfampos)
	if pfampos not in dic_mutations: dic_mutations[name] = Mutation(acc, mutation, pfampos, mut_type)
	# if mut_type not in dic_hmm[pfampos].mut_types:
	# 	dic_hmm[pfampos].mut_types[mut_type] = []
	# dic_hmm[pfampos].mut_types[mut_type].append(mutation)


'''Load ptm information for a pfam_pos'''
mycursor.execute("SELECT acc, uniprotpos, pfampos, ptmtype FROM ptms\
		 		where pfampos!='-'")
ptms = mycursor.fetchall()
dic_ptm_pfam = {}
dic_ptm = {}
dic_hmm = {}
for ptm_row in ptms:
	acc = ptm_row[0]
	uniprotpos = int(ptm_row[1])
	pfampos = int(ptm_row[2])
	ptm_type = ptm_row[3]
	for name in dic_mutations:
		position = int(dic_mutations[name].position)
		wtAA = dic_mutations[name].mutation[0]
		mutAA = dic_mutations[name].mutation[-1]
		if position == uniprotpos and acc == dic_mutations[name].acc:
			dic_mutations[name].ptm_types.append(ptm_type)
		elif uniprotpos in [position-1, position+1] and acc == dic_mutations[name].acc\
			and ptm_type in ['p', 'ac']:
			if wtAA in ['S', 'T', 'Y'] and mutAA in ['D', 'E']:
				dic_mutations[name].neighbour_pac_sites = 1
			elif wtAA in ['K'] and mutAA in ['Q']:
				dic_mutations[name].neighbour_pac_sites = 1
	# if pfampos not in dic_hmm: dic_hmm[pfampos] = HMM(pfampos)
	# dic_hmm[pfampos].ptm_types += 1

'''Load orthologs information'''
mycursor.execute("SELECT * FROM orth")
orth = mycursor.fetchall()
for orth_row in tqdm(orth):
	acc, wtAA, position = orth_row[0], orth_row[1], orth_row[2]
	position = int(position)
	num = 3
	for count, mutAA in enumerate(['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
								'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']):
		if mutAA != wtAA: continue
		orthscore = float(orth_row[num+count])
		for name in dic_mutations:
			if acc == dic_mutations[name].acc and position == dic_mutations[name].position:
				dic_mutations[name].orthscore = orthscore


'''Load HMM information for a pfam_pos'''
mycursor.execute("SELECT * FROM hmm")
hmm = mycursor.fetchall()
dic_hmm = {}
for hmm_row in hmm:
	pfampos, pfamaa = hmm_row[0], hmm_row[1]
	if pfampos == '-': continue
	pfampos = int(pfampos)
	num = 4
	for count, aa in enumerate(['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
								'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']):
		if aa == pfamaa.upper():
			bitscore = float(hmm_row[num+count])
			dic_hmm[pfampos] = bitscore
			break

# print (dic_hmm)
for name in dic_mutations:
	pfampos = dic_mutations[name].pfampos
	# print (pfampos)
	if pfampos not in dic_hmm: continue
	dic_mutations[name].bitscore = dic_hmm[pfampos]
	# print (name, dic_mutations[name].bitscore)

pfam_positions = [pfam_pos for pfam_pos in list(set(dic_hmm.keys()))]
pfam_positions.sort()

fig, ax = plt.subplots()

mut_types = ['constitutive-activation', 'increase', 'resistance', 'decrease', 'loss', 'neutral']
mut_types_colors = ['green', 'lightgreen', 'blue', 'lightcoral', 'red', '#F2E34C']
bottom = np.zeros(4)
for mut_type, mut_types_color in zip(mut_types, mut_types_colors):
	ptmsite = []
	for name in dic_mutations:
		if dic_mutations[name].mut_type != mut_type: continue
		mutation = dic_mutations[name].mutation
		num_ptms = len(dic_mutations[name].ptm_types)
		if num_ptms > 0:
			ptmsite.append(name)
	
	neighbour_pac_sites = []
	# for name in dic_mutations:
	# 	if dic_mutations[name].mut_type != mut_type: continue
	# 	mutation = dic_mutations[name].mutation
	# 	neighbour_pac_site = dic_mutations[name].neighbour_pac_sites
	# 	if neighbour_pac_site == 1:
	# 		neighbour_pac_sites.append(name)

	conserved_kinases = []
	conserved_orthologs = []
	for name in dic_mutations:
		if dic_mutations[name].mut_type != mut_type: continue
		if name in ptmsite: continue
		if name in neighbour_pac_sites: continue
		mutation = dic_mutations[name].mutation
		bitscore = dic_mutations[name].bitscore
		orthscore = dic_mutations[name].orthscore
		if bitscore <= 0.3:
			conserved_kinases.append(name)
		else:
			if orthscore >= 3.0:
				# print (name, orthscore, mut_type, dic_mutations[name].ptm_types, bitscore)
				conserved_orthologs.append(name)
		# if mut_type == 'constitutive-activation' and (orthscore >= 3.0 or bitscore <= 0.5):
		# 	print (name, orthscore, bitscore)
		if mut_type == 'resistance' and (orthscore >= 3.0 or bitscore <= 0.5):
			print (name, orthscore, bitscore)
	
	others = []
	for name in dic_mutations:
		if dic_mutations[name].mut_type != mut_type: continue
		if name in ptmsite: continue
		if name in neighbour_pac_sites: continue
		if name in conserved_kinases: continue
		if name in conserved_orthologs: continue
		mutation = dic_mutations[name].mutation
		others.append(name)

	row = [len(conserved_kinases), len(conserved_orthologs), len(ptmsite), len(others)]
	row = np.array(row)
	print (mut_type, row)
	width = 0.5
	p = ax.bar(['Conserved\nKinases', 'Conserved\nOrthologs', 'PTM', 'Others'], row, width,\
	    		label=mut_type, bottom=bottom, color=mut_types_color)
	bottom += row

ax.set_title("Distribution of mutations")
ax.legend(loc="upper center")
# plt.show()
plt.grid(linewidth=0.5, color='gray', linestyle='--')
plt.savefig('distribution_of_mutations.png', dpi=1000)
plt.savefig('distribution_of_mutations.eps', dpi=1000)
plt.savefig('distribution_of_mutations.svg', dpi=1000)