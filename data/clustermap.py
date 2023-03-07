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

## Class to store information
## about kinases
class Kinase:
	## initiate Object
	def __init__(self, acc):
		self.acc = acc
		self.pfam2seq = {}
		self.seq2pfam = {}
		self.fastaToAln = {}
		self.alnToFasta = {}
		self.mut_types = {'A':[], 'D':[], 'R':[]}

'''Load p-site information for a pfam_pos'''
dic_ptm_pfam = {}
for line in open('Kinase_psites4.tsv', 'r'):
	if line[0] == '#': continue
	if line.split('\t')[2] != 'Pkinase': continue
	ptm_type = line.split('\t')[3].split('-')[1]
	pfam_pos = line.split('\t')[4]
	if ptm_type == 'p':
		if pfam_pos not in dic_ptm_pfam:
			dic_ptm_pfam[pfam_pos] = 0
		else:
			dic_ptm_pfam[pfam_pos] += 1

kinases = {}
'''Store Uniprot accessions mappings'''
for line in gzip.open('../data/humanKinasesHmmsearchMappings2.tsv.gz', 'rt'):
	if line[0] == '#': continue
	acc = line.split('\t')[0].split('|')[1]
	uniprot_pos = line.split('\t')[2]
	pfam_pos = line.split('\t')[4].replace('\n', '')
	if acc not in kinases: kinases[acc] = Kinase(acc)
	kinases[acc].seq2pfam[uniprot_pos] = pfam_pos
	kinases[acc].pfam2seq[pfam_pos] = uniprot_pos

'''Load resistant mutation data'''
for line in gzip.open('../KA/resistant_mutations_Mar_2023.tsv.gz', 'rt'):
	if line[0] == '#': continue
	uniprot_acc = line.split('\t')[2]
	cosmic_pos = line.split('\t')[4]
	uniprot_pos = line.split('\t')[5].replace('\n', '')
	cosmic_mutation = line.split('\t')[1]
	uniprot_mutation = cosmic_mutation[0] + uniprot_pos + cosmic_mutation[-1]
	if uniprot_mutation not in kinases[uniprot_acc].mut_types['R']:
		kinases[uniprot_acc].mut_types['R'].append(uniprot_mutation)

print (kinases['P36507'].mut_types['R'])

'''Load act/dact mutation data'''
for line in open('../AK_mut_w_sc_feb2023/act_deact_v2.tsv', 'r'):
	if line.split()[0] == 'uniprot_name': continue
	uniprot_acc = line.split('\t')[1]
	uniprot_pos = line.split('\t')[3]
	wtAA = line.split('\t')[2].replace(',', '')
	mutAA = line.split('\t')[4].replace(',', '')
	if len(wtAA) > 1 or len(mutAA) > 1:
		continue
	uniprot_mutation = wtAA + uniprot_pos + mutAA
	mut_type = line.split('\t')[5]
	if mut_type == 'A':
		kinases[uniprot_acc].mut_types['A'].append(uniprot_mutation)
	else:
	    kinases[uniprot_acc].mut_types['D'].append(uniprot_mutation)

print (kinases['P36507'].mut_types['R'])

col_colors = []
for pfam_pos in range(1,265):
	if str(pfam_pos) not in dic_ptm_pfam:
		col_colors.append('white')
		continue
	if dic_ptm_pfam[str(pfam_pos)] >= 1:
		col_colors.append('grey')
	else:
		col_colors.append('white')

data = []
for mut_type in ['A', 'D', 'R']:
	row = []
	for pfam_pos in range(1,265):
		value = 0
		for acc in kinases:
			for mutation in kinases[acc].mut_types[mut_type]:
				# print (mutation	)
				uniprot_position = mutation[1:-1]
				wtAA = mutation[0]
				mutAA = mutation[-1]
				hmm_position = kinases[acc].seq2pfam[uniprot_position]
				if str(hmm_position) == str(pfam_pos):
					value += 1
					break
			if value  == 1:
				break
		row.append(value)
	data.append(row)

# print (data)

df = pd.DataFrame(
				data,
				index = ['A', 'D', 'R'],
				columns = [str(pfam_pos) for pfam_pos in range(1,265)])
df = df.loc[:, (df != 0).any(axis=0)]

print (df)

palette = sns.light_palette("seagreen", as_cmap=True)
g = sns.clustermap(df, cmap=palette, col_cluster=False, xticklabels=True, col_colors=col_colors)
#plt.savefig('clusterMap.jpeg', format="jpeg", dpi=500)
plt.show()
sys.exit()
## Store kinase information in
## the dictionary kinases
kinases = {}; allAlnPos = []
for line in open('../KA/hmmAlignmentMappings3.tsv', 'r'):
	if line[0] != '#':
		name = line.split('\t')[0]
		if name not in kinases:
			kinases[name] = Kinases(name)
		alnPos, fastaPos, fastaAA = int(line.split('\t')[1]), int(line.split('\t')[2]), str(line.replace('\n', '').split('\t')[3])
		if fastaAA.isupper():
			allAlnPos.append(alnPos)
			kinases[name].fastaToAln[fastaPos] = alnPos
			kinases[name].alnToFasta[alnPos] = fastaPos

def find_name(accGiven):
	for kinase in kinases:
		if accGiven in kinase:
			for name in kinase.split('_'):
				#print(name, accGiven)
				if name == accGiven:
					return kinase
	
	return None

## activating/inactivating mutations from UniProt
for line in open('../KA/act_deact_mut_for_scores_fin.tsv', 'r'):
	if line.split('\t')[0] != 'uniprot_id':
		acc = line.split('\t')[0]
		name = find_name(acc)
		mutation = line.split('\t')[1] + line.split('\t')[2] + line.split('\t')[3]
		status = line.replace('\n', '').split('\t')[-1]
		if status == 'A':
			kinases[name].act.append(mutation)
		else:
			kinases[name].deact.append(mutation)

## Resistance mutations from COSMIC
for line in open('../KA/resistance_mutations_w_scores_aligned_fin.tsv', 'r'):
	if line.split('\t')[0] != 'Gene.Name':
		name = line.split('\t')[0]
		mutation = line.split('\t')[1]
		kinases[name].resistance.append(mutation)

def findActDeactRes(alnPos):
	'''
	function to extract number of act, deact and resistance there are for
	the given alnPos
	'''
	act, deact, res = 0, 0, 0
	for kinase in kinases:
		for num, dic in enumerate([kinases[kinase].act, kinases[kinase].deact, kinases[kinase].resistance]):
			for mutation in dic:
				if mutation[-3:] == 'del':
					continue
				# print (mutation[-3:], alnPos, kinase, num)
				fastaPos = int(mutation[1:-1])
				if fastaPos not in kinases[kinase].fastaToAln:
					continue
				if alnPos == kinases[kinase].fastaToAln[fastaPos]:
					if num == 0:
						act += 1
					elif num == 1:
						deact += 1
					elif num == 2:
						res += 1
	return (act, deact, res)

allAlnPos = list(set(allAlnPos))
selected_allAlnPos = []
data = []
for count, alnPos in enumerate(allAlnPos):
	act, deact, res = findActDeactRes(alnPos)
	#print (count+1, len(allAlnPos))
	actBinary = 1 if act > 0 else 0
	deactBinary = 1 if deact > 0 else 0
	resBinary = 1 if res > 0 else 0
	if actBinary == 1 or deactBinary == 1 or resBinary == 1:
		selected_allAlnPos.append(alnPos)
		row=[]
		row.append(actBinary)
		row.append(deactBinary)
		row.append(resBinary)
		data.append(row)

#print (data)
df = pd.DataFrame(data, columns=['Activating', 'Deactivating', 'Resistance'])
print (df)
df.index = selected_allAlnPos
palette = sns.light_palette("seagreen", as_cmap=True)
g = sns.clustermap(df, cmap=palette)
#plt.savefig('clusterMap.jpeg', format="jpeg", dpi=500)
plt.show()
sys.exit()
