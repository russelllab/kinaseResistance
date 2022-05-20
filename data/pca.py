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

## Sequence features

homologyPath = '/net/home.isilon/ds-russell/mechismoX/analysis/alignments/data/HUMAN/orthologs_only/'

## Orthologs all
for acc in accs:
	orthoFile = homologyPath + acc[:4] + '/' + acc + '_orth.scores.txt.gz'
	for line in gzip.open(orthoFile, 'rt'):
		mut = line.split()[0].split('/')[1]
		score = float(line.split()[4])
		check = line.split()[1]
		#print (mut, diff)
		if mut in accs[acc].muts and check != '?':
			accs[acc].muts[mut].ortho = score

## Paralogs exlclusive
for acc in accs:
	orthoFile = homologyPath + acc[:4] + '/' + acc + '_excl_para.scores.txt.gz'
	for line in gzip.open(orthoFile, 'rt'):
		mut = line.split()[0].split('/')[1]
		score = float(line.split()[4])
		check = line.split()[1]
		#print (mut, diff)
		if mut in accs[acc].muts and check != '?':
			accs[acc].muts[mut].exclPara = score

## Paralogs specific
for acc in accs:
	orthoFile = homologyPath + acc[:4] + '/' + acc + '_spec_para.scores.txt.gz'
	for line in gzip.open(orthoFile, 'rt'):
		mut = line.split()[0].split('/')[1]
		score = float(line.split()[4])
		check = line.split()[1]
		#print (mut, diff)
		if mut in accs[acc].muts and check != '?':
			accs[acc].muts[mut].specPara = score

## Homologs
for acc in accs:
	orthoFile = homologyPath + acc[:4] + '/' + acc + '_all_homs.scores.txt.gz'
	for line in gzip.open(orthoFile, 'rt'):
		mut = line.split()[0].split('/')[1]
		score = float(line.split()[4])
		#if acc == 'P35968':
		#	print (mut, score)
		check = line.split()[1]
		#print (mut, diff)
		if mut in accs[acc].muts and check != '?':
			accs[acc].muts[mut].homo = score

## Structure features
AFpath = '/net/home.isilon/ds-russell/alphafold/HUMAN/'

for acc in accs:
	StrucFile = AFpath + '/AF-' + acc + '-F1-model_v1.dssp-scores.gz'
	for line in gzip.open(StrucFile, 'rt'):
		wt = line.split()[2]
		sub = line.split()[10]
		#print (line.split()[0])
		pos = int(line.split()[0])
		mut = wt + str(pos) + sub
		torsion = float(line.split()[18])
		sec = float(line.split()[22])
		bur = float(line.split()[26])
		accessibility = float(line.split()[30])
		#if acc == 'Q9UM73':
		#	print (wt, pos, sub, torsion, sec, bur, accessibility)
		if mut in accs[acc].muts:
			accs[acc].muts[mut].torsion = torsion
			accs[acc].muts[mut].burial = bur
			accs[acc].muts[mut].secondary = sec
			accs[acc].muts[mut].accessibility = accessibility

for acc in accs:
	StrucFile = AFpath + '/AF-' + acc + '-F1-model_v1.iupred.gz'
	for line in gzip.open(StrucFile, 'rt'):
		wt = line.split()[2]
		sub = line.split()[3]
		#print (line.split()[0])
		pos = int(line.split()[0])
		mut = wt + str(pos) + sub
		iupred = float(line.replace('\n', '').split()[-1])
		if acc == 'Q9UM73':
			print (wt, pos, sub, iupred)
		if mut in accs[acc].muts:
			accs[acc].muts[mut].iupred = iupred

for acc in accs:
	StrucFile = AFpath + '/AF-' + acc + '-F1-model_v1.mech_intra.gz'
	for line in gzip.open(StrucFile, 'rt'):
		if line.split()[0] == 'MECH':
			wt = line.split()[2]
			sub = line.split()[3]
			#print (line.split()[0])
			pos = int(line.split()[1])
			mut = wt + str(pos) + sub
			mechIntra = float(line.replace('\n', '').split()[6])
			if acc == 'Q9UM73':
				print (wt, pos, sub, mechIntra)
			if mut in accs[acc].muts:
				accs[acc].muts[mut].mechIntra = mechIntra
		'''
		else:
			wt = line.split()[2]
			sub = line.split()[3]
			#print (line.split()[0])
			pos = int(line.split()[1])
			mut = wt + str(pos) + sub
			ncont = float(line.replace('\n', '').split()[6])
			nres = float(line.replace('\n', '').split()[6])
			if acc == 'Q9UM73':
				print (wt, pos, sub, nres, ncont)
			if mut in accs[acc].muts:
				accs[acc].muts[mut].ncont = ncont
				accs[acc].muts[mut].nres = nres
		'''
									
#print (accs['P35968'].muts['Y1054F'].homo)
#sys.exit()
data = []
for acc in accs:
	for mut in accs[acc].muts:
		row = []
		if accs[acc].muts[mut].ortho!=None and accs[acc].muts[mut].specPara!=None and accs[acc].muts[mut].exclPara!=None:
			neutral = accs[acc].muts[mut].neutral
			#print(accs[acc].muts[mut].ortho)
			## Sequence
			ortho = float(accs[acc].muts[mut].ortho)
			specPara = float(accs[acc].muts[mut].specPara)
			exclPara = float(accs[acc].muts[mut].exclPara)
			homo = float(accs[acc].muts[mut].homo)
			
			## Structure
			torsion = float(accs[acc].muts[mut].torsion)
			secondary = float(accs[acc].muts[mut].secondary)
			burial = float(accs[acc].muts[mut].burial)
			iupred = float(accs[acc].muts[mut].iupred)
			accessibility = float(accs[acc].muts[mut].accessibility)
			mechIntra = float(accs[acc].muts[mut].mechIntra)
			#print (homo)
			
			row = []
			row.append(acc)
			row.append(mut)
			if accs[acc].muts[mut].resistance == True:
				row.append('Resistance')
				row.append('blue')
			elif accs[acc].muts[mut].activating == True:
				row.append('Activating')
				row.append('green')
			elif accs[acc].muts[mut].deactivating == True:
				row.append('Deactivating')
				row.append('red')
			elif accs[acc].muts[mut].neutral >= 20:
				#continue
				row.append('Neutral')
				row.append('grey')
			else:
				continue
				row.append(None)
			row.append(neutral)
			
			## Sequence
			row.append(ortho)
			row.append(exclPara)
			row.append(specPara)
			row.append(homo)
			
			## Structure
			row.append(torsion)
			row.append(secondary)
			row.append(burial)
			row.append(accessibility)
			row.append(iupred)
			row.append(mechIntra)
			
			data.append(row)
			#print (row)
#sys.exit()
df = pd.DataFrame(data, columns=['Accs',
								'Mutations',
								'Type',
								'Color',
								'HomoSamples',
								'Ortho',
								'ExclPara',
								'SpecPara',
								'AllHomo',
								'Torsion',
								'Secondary',
								'Burial',
								'Accessibility',
								'IUPred',
								'MechIntra'
								])
#df = pd.DataFrame(data, columns=['Mutations', 'Label', 'Type', 'Score'])
#newDF = df[df['Label'] == 'neutral']
X = df[['Ortho',
		'ExclPara',
		'SpecPara',
		'AllHomo',
		'Torsion',
		'Secondary',
		'Burial',
		'Accessibility',
		'IUPred',
		'MechIntra'
		]].values
		
pca = PCA(n_components=2)
Xt = pca.fit_transform(X)
print(pca.explained_variance_ratio_)
print(pca.singular_values_)
print(pca.components_)

#Xt = TSNE(n_components=2, learning_rate='auto', init='random').fit_transform(X)
dfNew = pd.concat([df['Type'], df['Color'], pd.DataFrame(Xt, columns = ['PC1', 'PC2'])], axis=1)
print (dfNew)
'''
ax = plt.axes(projection='3d')
zdata = dfNew['PC1'].values
xdata = dfNew['PC2'].values
ydata = dfNew['PC3'].values
cdata = dfNew['Color'].values
ax.scatter3D(xdata, ydata, zdata, c=cdata, cmap='Greens');
'''

ax = sns.scatterplot(data=dfNew, x="PC1", y="PC2", hue="Type",
				palette={'Activating': 'green', 'Deactivating':'red', 'Resistance':'blue', 'Neutral':'grey'},
				hue_order=['Activating', 'Deactivating', 'Resistance', 'Neutral'],
				)

#g = sns.FacetGrid(df, col="Type", hue="Label")
#g.map_dataframe(sns.violinplot, x="Label", y="Score")
#sns.scatterplot(data=df, x="Ortho", y="SpecPara", hue="Label", palette="deep")
'''
g = sns.catplot(x="Label", y="Score",
                col="Type",
                data=df, kind="box",
                height=4, aspect=.7);
'''
#ax = sns.boxplot(x="Label", y="Ortho", data=df)
ax.set_xlabel('PC1', fontsize=10)
ax.set_xlabel('PC2', fontsize=10)
#plt.savefig('pcaPlot.jpeg')
plt.show()
