#!/usr/bin/env python3
# coding: utf-8

import os, sys, gzip
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

dic = {}
for line in open('../KA/hmmAlignment.aln', 'r'):
	if line.split()!=[]:
		if line[0] != '#' and line.split()[0] != '//' and line.split()[0]!='CLUSTAL':
			name = line.split()[0]
			if 'UniProt' in name:
				acc = name.split('_UniProt')[0].split('_')[-1]
				seq = line.split()[1]
				#print (acc, seq)
				if acc not in dic:
					dic[acc] = ''
				dic[acc] += seq

#print (dic['FGFR2_P21802_UniProt'][:1080])
#print (dic['P21802'][:1080])
#sys.exit()
l = '#Name\tAlnPosition\tFastaPosition\tFastaAA\n'
fasta2Alignment = {}
for name in dic:
	aln = 0; fasta = 0
	if name not in fasta2Alignment:
		fasta2Alignment[name] = {}
	#print (name, dic[name])
	for char in dic[name]:
		#print (name + '\t' + str(aln) + '\t' + str(fasta) + '\t' + str(char) + '\t' + str(len(dic[name].replace('.', '').replace('-',''))) + '\n')
		if char!='.' and char!='-':
			l += str(name) + '\t' + str(aln) + '\t' + str(fasta) + '\t' + str(char) + '\n'
			fasta2Alignment[name][fasta] = aln
			fasta += 1
		aln += 1


class accessions:
    def __init__(self, acc, gene):
        self.acc = acc
        self.gene = gene
        self.muts = {}

class mutations(accessions):
	def __init__(self, mut):
		self.mut = mut
		self.ortho = None
		self.exclPara = None
		self.specPara = None
		self.homo = None
		self.activating = None
		self.deactivating = None
		self.neutral = None
		self.resistance = None
		self.torsion = None
		self.secondary = None
		self.accessibility = None
		self.burial = None
		self.iupred = None
		self.mechIntra = None

class alignments(accessions):
	def __init__(self, alnPos):
		self.alnPos = alnPos
		self.activating = None
		self.deactivating = None
		self.neutral = None
		self.resistance = None

accs = {}; aln = {}
## Activating mutations from UniProt
for line in open('../KA/kinase_activating_mutations_uniprot_with_gene_name.csv', 'r'):
	if line.split(',')[0] != 'Gene':
		acc = line.split(',')[1]
		gene = line.split(',')[0]
		try:
			pos = int(line.split(',')[3])
			alnPos = fasta2Alignment[acc][pos]
			if acc not in accs:
				accs[acc] = accessions(acc, gene)
			mut = line.split(',')[2]+str(pos)+line.split(',')[4]

			if mut not in accs[acc].muts:
				accs[acc].muts[mut] = mutations(mut)
			accs[acc].muts[mut].activating = True

			if alnPos not in aln:
				aln[alnPos] = alignments(alnPos)
			aln[alnPos].activating = True
		except:
			continue

## Deactivating mutations from UniProt
for line in open('../KA/kinase_deactivating_mutations_uniprot_with_gene_name.csv', 'r'):
	if line.split(',')[0] != 'Gene':
		acc = line.split(',')[1]
		gene = line.split(',')[0]
		pos = int(line.split(',')[3])
		alnPos = fasta2Alignment[acc][pos]
		if acc not in accs:
			accs[acc] = accessions(acc, gene)
		mut = line.split(',')[2]+str(pos)+line.split(',')[4]

		if mut not in accs[acc].muts:
			accs[acc].muts[mut] = mutations(mut)
		accs[acc].muts[mut].deactivating = True

		if alnPos not in aln:
			aln[alnPos] = alignments(alnPos)
		aln[alnPos].deactivating = True

## gnomAD data
gnomadPath = '/net/home.isilon/ds-russell/mechismoX/analysis/mutations/data/VLatest/'
for acc in accs:
	if os.path.isfile(gnomadPath + acc + '.annotations.txt.gz'):
		for line in gzip.open(gnomadPath + acc + '.annotations.txt.gz', 'rt'):
			#print (line.split('\t')[0].split())
			if line.split('\t')[0].split() != ['UniProtID']:
				hom = line.split('\t')[15].rstrip()
				if hom != '-':
					hom = int(hom)
					if hom > 0:
						pos = int(line.split('\t')[2].rstrip())
						alnPos = fasta2Alignment[acc][pos]
						mut = line.split('\t')[1].rstrip()+str(pos)+line.split('\t')[3].rstrip()
						#print (acc, accs[acc].gene, int(hom), mut)
						if mut not in accs[acc].muts:
							accs[acc].muts[mut] = mutations(mut)
						accs[acc].muts[mut].neutral = hom

						if alnPos not in aln:
							aln[alnPos] = alignments(alnPos)
						aln[alnPos].neutral = hom

	else:
		print (acc, accs[acc].gene)

## Resistance mutations from COSMIC
for line in open('../KA/kinase_resistant_mutation_count_by_ID_sample.csv', 'r'):
	if line.split(',')[0] != 'Gene':
		gene = line.split(',')[0]
		for acc in accs:
			if accs[acc].gene == gene:
				pos = int(line.split(',')[2])
				alnPos = fasta2Alignment[acc][pos]
				mut = line.split(',')[1]+str(pos)+line.split(',')[3]
				if mut not in accs[acc].muts:
					accs[acc].muts[mut] = mutations(mut)
				accs[acc].muts[mut].resistance = True

				if alnPos not in aln:
					aln[alnPos] = alignments(alnPos)
				aln[alnPos].resistance = True
				break

data = []
for alnPos in aln:
	print (alnPos, aln[alnPos].activating, aln[alnPos].deactivating, aln[alnPos].resistance)
	row = []
	if aln[alnPos].activating == True:
		row.append(1)
	else:
		row.append(0)
	if aln[alnPos].deactivating == True:
		row.append(1)
	else:
		row.append(0)
	if aln[alnPos].resistance == True:
		row.append(1)
	else:
		row.append(0)
	data.append(row)

print (data)
df = pd.DataFrame(data, columns=['Activating', 'Deactivating', 'Resistance'])
palette = sns.light_palette("seagreen", as_cmap=True)
g = sns.clustermap(df, cmap=palette)
plt.savefig('clusterMap.jpeg', format="jpeg", dpi=1000)
plt.show()
sys.exit()
