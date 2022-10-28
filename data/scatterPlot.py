#!/usr/bin/env python3
# coding: utf-8

import os, sys, gzip
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

dic = {}
for line in open('../KA/hmmAlignment.aln', 'r'):
    if line.split()!=[]:
        if line[0] != '#' and line.split()[0] != '//':
            name = line.split()[0]
            seq = line.split()[1]
            #print (name, seq)
            if name not in dic:
                dic[name] = ''
            dic[name] += seq

#print (dic['FGFR2_P21802_UniProt'][:1080])
#sys.exit()
l = '#Name\tAlnPosition\tFastaPosition\tFastaAA\n'
fasta2Alignment = {}
for name in dic:
	aln = 1; fasta = 1
	if name not in fasta2Alignment:
		fasta2Alignment[name] = {}
	#print (name, dic[name])
	for char in dic[name]:
		#print (name + '\t' + str(aln) + '\t' + str(fasta) + '\t' + str(char) + '\t' + str(len(dic[name].replace('.', '').replace('-',''))) + '\n')
		if char!='.' and char!='-':
			l += str(name) + '\t' + str(aln) + '\t' + str(fasta) + '\t' + str(char) + '\n'
			fasta += 1
			fasta2Alignment[name][fasta] = aln
		aln += 1
        

class accessions:
    def __init__(self, acc, gene):
        self.acc = acc
        self.gene = gene
        self.muts = {}

class mutations():
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
        
accs = {}
## Activating mutations from UniProt
for line in open('../KA/kinase_activating_mutations_uniprot_with_gene_name.csv', 'r'):
    if line.split(',')[0] != 'Gene':
        acc = line.split(',')[1]
        gene = line.split(',')[0]
        if acc not in accs:
            accs[acc] = accessions(acc, gene)
        mut = line.split(',')[2]+line.split(',')[3]+line.split(',')[4]
        if mut not in accs[acc].muts:
            accs[acc].muts[mut] = mutations(mut)
            accs[acc].muts[mut].activating = True

## Deactivating mutations from UniProt
for line in open('../KA/kinase_deactivating_mutations_uniprot_with_gene_name.csv', 'r'):
    if line.split(',')[0] != 'Gene':
        acc = line.split(',')[1]
        gene = line.split(',')[0]
        if acc not in accs:
            accs[acc] = accessions(acc, gene)
        mut = line.split(',')[2]+line.split(',')[3]+line.split(',')[4]
        if mut not in accs[acc].muts:
            accs[acc].muts[mut] = mutations(mut)
            accs[acc].muts[mut].deactivating = True

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
						mut = line.split('\t')[1].rstrip()+line.split('\t')[2].rstrip()+line.split('\t')[3].rstrip()
						#print (acc, accs[acc].gene, int(hom), mut)
						if mut not in accs[acc].muts:
							accs[acc].muts[mut] = mutations(mut)
						
						accs[acc].muts[mut].neutral = hom
						
	else:
		print (acc, accs[acc].gene)

## Resistance mutations from COSMIC
for line in open('../KA/kinase_resistant_mutation_count_by_ID_sample.csv', 'r'):
	if line.split(',')[0] != 'Gene':
		gene = line.split(',')[0]
		for acc in accs:
			if accs[acc].gene == gene:
				mut = line.split(',')[1]+line.split(',')[2]+line.split(',')[3]
				if mut not in accs[acc].muts:
					accs[acc].muts[mut] = mutations(mut)
					accs[acc].muts[mut].resistance = True
				break

homologyPath = '/net/home.isilon/ds-russell/mechismoX/analysis/alignments/data/HUMAN/orthologs_only/'

## Orthologs all
for acc in accs:
	orthoFile = homologyPath + acc[:4] + '/' + acc + '_orth.scores.txt.gz'
	for line in gzip.open(orthoFile, 'rt'):
		mut = line.split()[0].split('/')[1]
		score = float(line.split()[4])
		#print (mut, diff)
		if mut in accs[acc].muts:
			accs[acc].muts[mut].ortho = score

## Paralogs exlclusive
for acc in accs:
	orthoFile = homologyPath + acc[:4] + '/' + acc + '_excl_para.scores.txt.gz'
	for line in gzip.open(orthoFile, 'rt'):
		mut = line.split()[0].split('/')[1]
		score = float(line.split()[4])
		#print (mut, diff)
		if mut in accs[acc].muts:
			accs[acc].muts[mut].exclPara = score

## Paralogs specific
for acc in accs:
	orthoFile = homologyPath + acc[:4] + '/' + acc + '_spec_para.scores.txt.gz'
	for line in gzip.open(orthoFile, 'rt'):
		mut = line.split()[0].split('/')[1]
		score = float(line.split()[4])
		#print (mut, diff)
		if mut in accs[acc].muts:
			accs[acc].muts[mut].specPara = score

## Homologs
for acc in accs:
	orthoFile = homologyPath + acc[:4] + '/' + acc + '_all_homs.scores.txt.gz'
	for line in gzip.open(orthoFile, 'rt'):
		mut = line.split()[0].split('/')[1]
		score = float(line.split()[4])
		#print (mut, diff)
		if mut in accs[acc].muts:
			accs[acc].muts[mut].homo = score

data = []
for acc in accs:
	for mut in accs[acc].muts:
		row = []
		neutral = accs[acc].muts[mut].neutral
		ortho = accs[acc].muts[mut].ortho
		specPara = accs[acc].muts[mut].specPara
		exclPara = accs[acc].muts[mut].exclPara
		homo = accs[acc].muts[mut].homo
		
		row = []
		row.append(mut)
		if accs[acc].muts[mut].resistance == True:
			row.append('resistance')
		elif accs[acc].muts[mut].activating == True:
			row.append('activating')
		elif accs[acc].muts[mut].deactivating == True:
			row.append('deactivating')
		elif accs[acc].muts[mut].neutral > 50:
			row.append('neutral')
		else:
			continue
			row.append(None)
		row.append(neutral)		
		row.append(ortho)
		row.append(exclPara)
		row.append(specPara)
		row.append(hom)
		data.append(row)
		'''
		for name, value in zip(['Orthologs', 'ExclParalogs', 'SpecParalogs', 'AllHomologs'],[ortho, exclPara, specPara, homo]):
			row = []
			row.append(mut)
			if accs[acc].muts[mut].resistance == True:
				row.append('resistance')
			elif accs[acc].muts[mut].activating == True:
				row.append('activating')
			elif accs[acc].muts[mut].deactivating == True:
				row.append('deactivating')
			elif accs[acc].muts[mut].neutral > 5:
				row.append('neutral')
			else:
				continue
				row.append(None)
			
			#row.append(neutral)		
			#row.append(ortho)
			#row.append(exclPara)
			#row.append(specPara)
			#row.append(hom)
			
			row.append(name)
			row.append(value)
			data.append(row)
		'''
		if ortho != None and neutral != None:
			print (acc, mut, accs[acc].muts[mut].neutral, accs[acc].muts[mut].ortho, accs[acc].muts[mut].activating)

df = pd.DataFrame(data, columns=['Mutations', 'Label', 'HomoSamples', 'Ortho', 'ExclPara', 'SpecPara', 'AllHomo'])
#df = pd.DataFrame(data, columns=['Mutations', 'Label', 'Type', 'Score'])
print (df)
#g = sns.FacetGrid(df, col="Type", hue="Label")
#g.map_dataframe(sns.violinplot, x="Label", y="Score")
sns.scatterplot(data=df, x="Ortho", y="SpecPara", hue="Label", palette="deep")
'''
g = sns.catplot(x="Label", y="Score",
                col="Type",
                data=df, kind="box",
                height=4, aspect=.7);
'''
#ax = sns.boxplot(x="Label", y="Ortho", data=df)
#plt.savefig('scatterplot.jpeg')
plt.show()
