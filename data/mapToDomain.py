import os, sys, gzip

for line in gzip.open('CosmicResistanceMutations.tsv.gz', 'rt'):
	if 'Sample Name' not in line:
		if 'Confirmed Somatic' in line:
			gene = line.split('\t')[2]
			mutation = line.split('\t')[9]
			print (gene, mutation)
		
