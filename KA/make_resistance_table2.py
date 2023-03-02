#!/usr/bin/env python3
# coding: utf-8

'''
Script to make resistance table
'''

import os, sys, gzip

COSMIC_FILE = '../AK_mut_w_sc_feb2023/res_mut_old.tsv'

class Genes:
    def __init__(self, gene) -> None:
        self.gene = gene
        self.mutations = {}
        self.mappings = {}
    
    def store_mappings(self, seq_pos, pfam_pos):
        self.mappings[seq_pos] = pfam_pos

class Mutations:
    def __init__(self, mutation, zygosity, census_gene, transcript, canonical_acc) -> None:
        self.mutation = mutation
        self.transcript = transcript
        self.cannonical_acc = cannonical_acc
        self.zygosity = [zygosity]
        self.census_gene = [census_gene]
        # self.resistant_to_inhib = [resistant_to_inhib]
        # self.sample_id = [sample_id]

def get_acc(gene, mutation_address):
    canonical_acc = mutation_address.cannonical_acc
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
    if acc == canonical_acc:
        print (acc, gene, mutation_address.mutation, mutation_address.transcript)
        # sys.exit()
    if '-' in acc:
        if acc.split('-')[1] == '1':
            acc = acc.split('-')[0]
    # print (acc)
    return acc

genes = {}
for line in open(COSMIC_FILE, 'r'):
    if line.split('\t')[0] == 'Gene.Name':
        ## Ignore the header
        continue
    # sample_id = line.split('\t')[1]
    gene = line.split('\t')[0]
    # gene = gene if '_' not in gene else gene.split('_')[0]
    transcript = line.split('\t')[15]
    census_gene = line.split('\t')[1]
    # resistant_to_inhib = line.split('\t')[5]
    zygosity = line.split('\t')[14]
    mutation = line.split('\t')[3]
    cannonical_acc = line.split('\t')[17]
    if '?' not in mutation:
    # if '?' not in mutation:
        if gene not in genes:
            genes[gene] = Genes(gene)
        if mutation not in genes[gene].mutations:
            genes[gene].mutations[mutation] = Mutations(mutation,
                                                        zygosity,
                                                        census_gene,
                                                        transcript,
                                                        cannonical_acc)
        else:
            genes[gene].mutations[mutation].zygosity.append(zygosity)
            genes[gene].mutations[mutation].census_gene.append(census_gene)
            # genes[gene].mutations[mutation].resistant_to_inhib.append(resistant_to_inhib)
            # genes[gene].mutations[mutation].sample_id.append(sample_id)

'''
Adding mapping information of cosmic genes
'''
for line in gzip.open('../data/cosmic_kinases_Mappings.tsv.gz', 'rt'):
    if line[0] == '#':
        continue
    gene = line.split('\t')[0].split(':')[0]
    seq_pos = line.split('\t')[2]
    pfam_pos = line.replace('\n', '').split('\t')[4]
    genes[gene].store_mappings(seq_pos, pfam_pos)

print (genes['EGFR_ENST00000454757'].mappings)

'''
Get mapping information of UniProt accs
'''
uniprot_mappings = {}
for line in gzip.open('../data/humanKinasesHmmsearchMappings2.tsv.gz', 'rt'):
    if line[0] == '#':
        continue
    acc = line.split('\t')[0].split('|')[1]
    seq_pos = line.split('\t')[2]
    pfam_pos = line.replace('\n', '').split('\t')[4]
    if acc not in uniprot_mappings:
        uniprot_mappings[acc] = {}
    else:
        uniprot_mappings[acc][seq_pos] = pfam_pos
## complete here 
sys.exit()
# print (genes['ABL1_ENST00000318560'].mutations['T315I'].cannonical_acc)
# sys.exit()
out_text = '#gene\tuniprot_acc\tmutation\tzygosity\tnum_samples\tcensus_gene\tresistant_to_inhib\n'
for gene in genes:
    # print (gene)
    for mutation in genes[gene].mutations:
        mutation_address = genes[gene].mutations[mutation]
        out_text += gene + '\t'
        # out_text += human_kinases[gene] + '\t'
        out_text += get_acc(gene, mutation_address) + '\n'
        # out_text += mutation.split('.')[1] + '\t'
        # out_text += ';'.join(list(set(mutation_address.zygosity))) + '\t'
        # out_text += str(len(list(set(mutation_address.sample_id)))) + '\t'
        # out_text += ';'.join(list(set(mutation_address.census_gene))) + '\t'
        # out_text += ';'.join(list(set(mutation_address.resistant_to_inhib))) + '\n'

# print (out_text)
sys.exit()
gzip.open('resistant_mutations_Nov22new.tsv.gz', 'wt').write(out_text)