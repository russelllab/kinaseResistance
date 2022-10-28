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
class Kinases:
	## initiate Object
	def __init__(self, name):
		self.name = name
		self.fastaToAln = {}
		self.alnToFasta = {}
		self.act = []
		self.deact = []
		self.resistance = []

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
