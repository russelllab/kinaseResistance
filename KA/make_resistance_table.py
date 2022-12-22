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
    def __init__(self, mutation, zygosity, census_gene, resistant_to_inhib, sample_id, transcript) -> None:
        self.mutation = mutation
        self.transcript = transcript
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

def get_acc(gene, mutation_address):
    canonical_acc = human_kinases[gene]
    acc = canonical_acc
    for line in gzip.open('UniProtFasta2/'+canonical_acc+'.txt.gz', 'rt'):
        if line[:2] != 'DR':
            continue
        if line.split()[1] != 'Ensembl;':
            continue
        if '[' not in line.split()[-1]:
            continue
        acc_found = line.split()[-1].replace('[','').replace(']', '')
        transcript_found = line.split(';')[1].replace(' ', '')
        # print (transcript_found)
        if mutation_address.transcript.split('.')[0] == transcript_found.split('.')[0]:
            acc = acc_found
            break
        # print (line.split()[-1])
        # print (line)
    # return canonical_acc
    if '-' in acc:
        if acc.split('-')[1] == '1':
            acc = acc.split('-')[0]
    # print (acc)
    return acc

genes = {}
for line in gzip.open(COSMIC_FILE, 'rt'):
    if line.split('\t')[:2] == ['Sample Name', 'Sample ID']:
        ## Ignore the header
        continue
    sample_id = line.split('\t')[1]
    gene = line.split('\t')[2]
    # gene = gene if '_' not in gene else gene.split('_')[0]
    transcript = line.split('\t')[3]
    census_gene = line.split('\t')[4]
    resistant_to_inhib = line.split('\t')[5]
    zygosity = line.split('\t')[21]
    mutation = line.split('\t')[9]
    if '?' not in mutation and gene in human_kinases:
    # if '?' not in mutation:
        if gene not in genes:
            genes[gene] = Genes(gene)
        if mutation not in genes[gene].mutations:
            genes[gene].mutations[mutation] = Mutations(mutation, zygosity, census_gene, resistant_to_inhib, sample_id, transcript)
        else:
            genes[gene].mutations[mutation].zygosity.append(zygosity)
            genes[gene].mutations[mutation].census_gene.append(census_gene)
            genes[gene].mutations[mutation].resistant_to_inhib.append(resistant_to_inhib)
            genes[gene].mutations[mutation].sample_id.append(sample_id)

# print (genes['MAP2K2'].mutations['p.V215E'].sample_id)
# sys.exit()
out_text = '#gene\tuniprot_acc\tmutation\tzygosity\tnum_samples\tcensus_gene\tresistant_to_inhib\n'
for gene in genes:
    # print (gene)
    for mutation in genes[gene].mutations:
        mutation_address = genes[gene].mutations[mutation]
        out_text += gene + '\t'
        # out_text += human_kinases[gene] + '\t'
        out_text += get_acc(gene, mutation_address) + '\t'
        out_text += mutation.split('.')[1] + '\t'
        out_text += ';'.join(list(set(mutation_address.zygosity))) + '\t'
        # out_text += ';'.join(list(set(mutation_address.sample_id))) + '\t'
        out_text += str(len(list(set(mutation_address.sample_id)))) + '\t'
        out_text += ';'.join(list(set(mutation_address.census_gene))) + '\t'
        out_text += ';'.join(list(set(mutation_address.resistant_to_inhib))) + '\n'

# print (out_text)
gzip.open('resistant_mutations_Nov22new.tsv.gz', 'wt').write(out_text)