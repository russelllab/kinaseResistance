#!/usr/bin/env python3

import os, sys, gzip

class Mutation:
    def __init__(self, mutation, mutation_type, gene, info, outcome):
        self.acc = mutation.split('/')[0]
        self.mut = mutation.split('/')[1]
        # self.mutation = mutation
        self.mutation_type = mutation_type
        self.gene = gene
        self.info = info
        self.outcome = outcome
    
    def show(self):
        # print(self.mutation, self.mutation_type, self.gene, self.info, self.outcome)
        return '\t'.join([
                        self.acc,
                        self.gene,
                        self.mut,
                        # self.mutation,
                        self.mutation_type,
                        self.info,
                        self.outcome
                        ])

dic = {}
for line in open('ana_set.tsv', 'r'):
    if line.startswith('UniProtAccMutation'): continue
    mutation= line.split('\t')[0]
    mutation_type = line.split('\t')[5]
    gene = line.split('\t')[3]
    info = line.split('\t')[9]
    outcome1 = line.split('\t')[10]
    outcome2 = line.split('\t')[13]
    if outcome2 == 'exclude': continue
    outcome = outcome2 if outcome2 != '' else outcome1
    if mutation not in dic:
        dic[mutation] = Mutation(mutation, mutation_type, gene, info, outcome)
        dic[mutation].show()

for line in open('missing_cases_annotated.tsv', 'r'):
    if line.startswith('UniProtAcc'): continue
    mutation= line.split('\t')[3]
    gene = line.split('\t')[1]
    info = line.split('\t')[6] + ' '+ line.split('\t')[7]
    mutation_type = line.split('\t')[5]
    outcome = line.split('\t')[9]
    if outcome == '': continue
    if mutation not in dic:
        dic[mutation] = Mutation(mutation, mutation_type, gene, info, outcome)
        dic[mutation].show()

for line in open('final_mined_RR_checked_checked-again.txt', 'r'):
    if line[0] == '#': continue
    mutation= line.split('\t')[0]
    mut = mutation.split('/')[1]
    # only missense mutations
    if mut[1:-1].isdigit() == False: continue
    if mut[0].isalpha() == False: continue
    if mut[-1].isalpha() == False: continue
    gene = line.split('\t')[1]
    info = line.split('\t')[3]+'; '+line.split('\t')[2]
    mutation_type = 'VARIANT'
    outcome = 'activating'
    if mutation not in dic:
        dic[mutation] = Mutation(mutation, mutation_type, gene, info, outcome)
        # dic[mutation].show()

text = 'UniProtAcc\tGene\tMutation\tMutationType\tInfo\tOutcome\n'
for mutation in dic:
    text += dic[mutation].show() + '\n'

gzip.open('ad_mutations.tsv.gz', 'wt').write(text)