#!/usr/bin/env python3
# coding: utf-8

'''
Script to make resistance table
'''

import os, sys, gzip

COSMIC_FILE = 'CosmicResistanceMutations_Nov22.tsv.gz'

class Genes:
    def __init__(self, gene) -> None:
        self.gene = gene
        self.mutations = {}

class Mutations:
    def __init__(self, mutation, zygosity, census_gene, resistant_to_inhib, sample_id) -> None:
        self.mutation = mutation
        self.zygosity = [zygosity]
        self.census_gene = [census_gene]
        self.resistant_to_inhib = [resistant_to_inhib]
        self.sample_id = [sample_id]

human_kinases = {}
for line in open('../data/humanKinases.fasta', 'r'):
    if line[0] != '>':
        continue
    gene = line.split('GN=')[1].split()[0]
    acc = line.split('|')[1]
    human_kinases[gene] = acc

genes = {}
for line in gzip.open(COSMIC_FILE, 'rt'):
    if line.split('\t')[:2] == ['Sample Name', 'Sample ID']:
        ## Ignore the header
        continue
    sample_id = line.split('\t')[1]
    gene = line.split('\t')[2]
    census_gene = line.split('\t')[4]
    resistant_to_inhib = line.split('\t')[5]
    zygosity = line.split('\t')[21]
    mutation = line.split('\t')[9]
    if '?' not in mutation and gene in human_kinases:
        if gene not in genes:
            genes[gene] = Genes(gene)
        if mutation not in genes[gene].mutations:
            genes[gene].mutations[mutation] = Mutations(mutation, zygosity, census_gene, resistant_to_inhib, sample_id)
        else:
            genes[gene].mutations[mutation].zygosity.append(zygosity)
            genes[gene].mutations[mutation].census_gene.append(census_gene)
            genes[gene].mutations[mutation].resistant_to_inhib.append(resistant_to_inhib)
            genes[gene].mutations[mutation].sample_id.append(sample_id)

# print (genes['MAP2K2'].mutations['p.V215E'].sample_id)
# sys.exit()
out_text = '#gene\tuniprot_acc\tmutation\tzygosity\tnum_samples\tcensus_gene\tresistant_to_inhib\n'
for gene in genes:
    for mutation in genes[gene].mutations:
        mutation_address = genes[gene].mutations[mutation]
        out_text += gene + '\t'
        out_text += human_kinases[gene] + '\t'
        out_text += mutation.split('.')[1] + '\t'
        out_text += ';'.join(list(set(mutation_address.zygosity))) + '\t'
        # out_text += ';'.join(list(set(mutation_address.sample_id))) + '\t'
        out_text += str(len(list(set(mutation_address.sample_id)))) + '\t'
        out_text += ';'.join(list(set(mutation_address.census_gene))) + '\t'
        out_text += ';'.join(list(set(mutation_address.resistant_to_inhib))) + '\n'

print (out_text)
open('resistant_mutations_Nov22.tsv.gz', 'wt').write(out_text)