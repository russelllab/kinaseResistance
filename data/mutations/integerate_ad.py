#!/usr/bin/env python3

import os, sys, gzip, re

class Mutation:
    def __init__(self, mutation, mutation_type, gene, info, pubmedIDs, outcome, source):
        self.acc = mutation.split('/')[0]
        self.mut = mutation.split('/')[1]
        # self.mutation = mutation
        self.mutation_type = mutation_type
        self.gene = gene
        self.info = info
        self.pubmedIDs = pubmedIDs
        self.outcome = outcome
        self.source = source
    
    def show(self):
        # print(self.mutation, self.mutation_type, self.gene, self.info, self.outcome)
        return '\t'.join([
                        self.acc,
                        self.gene,
                        self.mut,
                        # self.mutation,
                        self.mutation_type,
                        self.info,
                        self.pubmedIDs,
                        self.outcome,
                        self.source
                        ])

def extractInformation(text):
    '''
    extract pubmed information
    '''
    pubmed_pattern = r"PubMed:\s*(\d+)"
    pubmed_matches = re.findall(pubmed_pattern, text)
    pubmed = ','.join([match for match in pubmed_matches])
    pubmed = pubmed.rstrip('\n')

    return pubmed

info_dic = {}
for line in open('../tt3', 'r'):
    if line.startswith('\n'): continue
    if line.startswith('UniProt'): continue
    if line.startswith('Error:'): continue
    # print(line)
    mutation = line.split('\t')[3]
    info = line.split('\t')[5]+ ' ' + line.split('\t')[6]
    pubmedIDs = line.split('\t')[8].rstrip()
    if mutation not in info_dic:
        info_dic[mutation] = {}
        info_dic[mutation]['info'] = info
        info_dic[mutation]['pubmedIDs'] = pubmedIDs

dic = {}
for line in open('ana_set.tsv', 'r'):
    if line.startswith('UniProtAccMutation'): continue
    mutation= line.split('\t')[0]
    mutation_type = line.split('\t')[5]
    gene = line.split('\t')[3]
    # info = line.split('\t')[9]
    if mutation in info_dic:
        info = info_dic[mutation]['info']
        pubmedIDs = info_dic[mutation]['pubmedIDs']
    else:
        description = line.split('\t')[9].rstrip()
        info = description.split('"";/evidence')[0]
        pubmedIDs = extractInformation(description)
    outcome1 = line.split('\t')[10]
    outcome2 = line.split('\t')[13]
    if outcome2 == 'exclude': continue
    outcome = outcome2 if outcome2 != '' else outcome1
    # print (line)
    if mutation not in dic:
        dic[mutation] = Mutation(mutation, mutation_type, gene, info, pubmedIDs, outcome, 'UniProt')
        dic[mutation].show()

for line in open('missing_cases_annotated.tsv', 'r'):
    if line.startswith('UniProtAcc'): continue
    mutation= line.split('\t')[3]
    gene = line.split('\t')[1]
    info = line.split('\t')[5] + ' '+ line.split('\t')[6]
    pubmedIDs = line.split('\t')[8]
    mutation_type = line.split('\t')[5]
    outcome = line.split('\t')[9]
    if outcome == '': continue
    if mutation not in dic:
        dic[mutation] = Mutation(mutation, mutation_type, gene, info, pubmedIDs, outcome, 'UniProt')
        dic[mutation].show()

for line in gzip.open('final_mined_RR_checked_checked-again.txt.gz', 'rt'):
    if line[0] == '#': continue
    mutation= line.split('\t')[0]
    mut = mutation.split('/')[1]
    # only missense mutations
    if mut[1:-1].isdigit() == False: continue
    if mut[0].isalpha() == False: continue
    if mut[-1].isalpha() == False: continue
    gene = line.split('\t')[1]
    pubmedIDs = line.split('\t')[3]
    mutation_type = 'VARIANT'
    info = mutation_type + ' '+line.split('\t')[2]
    outcome = 'activating'
    if mutation not in dic:
        dic[mutation] = Mutation(mutation, mutation_type, gene, info, pubmedIDs, outcome, 'PubMed')
        # dic[mutation].show()
    else:
        dic[mutation].source += '+PubMed'

text = 'UniProtAcc\tGene\tMutation\tMutationType\tInfo\tPubMedID\tOutcome\tSource\n'
for mutation in dic:
    # print (mutation)
    text += dic[mutation].show() + '\n'
    # print (dic[mutation].show())

gzip.open('ad_mutations.tsv.gz', 'wt').write(text)