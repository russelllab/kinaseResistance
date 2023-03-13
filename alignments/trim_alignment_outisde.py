#!/usr/bin/env python3

'''
A script to take an alignment as input
and trim based on some paramerters
'''

import pandas as pd
import numpy as np
from tqdm import tqdm
import sys, gzip
import pickle

INPUT_FILE = 'humanKinases.clustal'
GAP_THRESHOLD = 0.9
START_ALN = 32178
END_ALN = 33550
WINDOW = 30

class Kinase:
    '''A class to store information about a kinase'''
    def __init__(self, acc, name, sequence):
        self.acc = acc
        self.name = name
        self.sequence = sequence
        self.fasta = ''
        self.aln2seq = {}
        self.seq2aln = {}
        self.mutations = {'A':[], 'D':[], 'R':[]}
    def convert2fasta(self) -> None:
        '''Convert the alignment to a fasta sequence'''
        self.fasta = self.sequence.replace('-', '').replace('.', '').upper()
    def make_aln2seq(self) -> None:
        '''Make a dictionary to convert alignment positions to sequence positions'''
        aln_pos = 0
        seq_pos = 0
        for aa in self.sequence:
            if aa in ['-', '.']:
                aln_pos += 1
            else:
                self.aln2seq[aln_pos] = seq_pos
                self.seq2aln[seq_pos] = aln_pos
                aln_pos += 1
                seq_pos += 1
    def find_fasta_position(self, aln_pos: int) -> str:
        '''Find the position in the fasta sequence'''
        while True:
            if aln_pos in self.aln2seq:
                return str(self.aln2seq[aln_pos] + 1)
            else:
                aln_pos += 1

kinases = {}
# open the file and save the sequences
for line in open(INPUT_FILE, 'r'):
    # ignore the header
    if line.startswith('CLUSTAL') or line.startswith('//'):
        continue
    # ignore the blank lines
    if line.strip() == '':
        continue
    # get the sequence and name
    sequence = line.split()[1]
    name = line.split()[0]
    acc = line.split()[0].split('|')[1]
    # add sequence to the object
    if acc not in kinases:
        kinases[acc] = Kinase(acc, name, sequence)
    else:
        kinases[acc].sequence += sequence

# convert the alignment to fasta
# and make the mapping dictionaries
aln_data = []
aln_accs = []
for num, acc in enumerate(kinases):
    # print (kinases['sp|Q92772|CDKL2_HUMAN'].fasta)
    # print (name)
    kinases[acc].convert2fasta()
    kinases[acc].make_aln2seq()
    aln_data.append(list(kinases[acc].sequence))
    aln_accs.append(acc)
    # if num == 2: break
    
df = pd.DataFrame(aln_data, index=aln_accs)
# print (df[range(START_ALN-WINDOW, END_ALN+WINDOW+1)])
df = df[range(START_ALN-WINDOW, END_ALN+WINDOW+1)]

trimmed_aln = 'CLUSTAL\n\n'
for acc, row in zip(aln_accs, df.to_numpy()):
    # print (name, ''.join(row))
    trimmed_aln += kinases[acc].name + '|'
    trimmed_aln += kinases[acc].find_fasta_position(START_ALN - WINDOW) + ' '
    trimmed_aln += ''.join(row) + '\n'

open('humanKinasesTrimmed.clustal', 'w').write(trimmed_aln)

'''Resistance Mutations'''
for line in gzip.open('../KA/resistant_mutations_Mar_2023.tsv.gz', 'rt'):
    if line[0] == '#':
        continue
    acc = line.split('\t')[2]
    cosmic_mutation = line.split('\t')[1]
    wtAA = cosmic_mutation[0]
    mutAA = cosmic_mutation[-1]
    uniprot_position = str(line.split('\t')[5].strip())
    uniprot_mutation = wtAA + uniprot_position + mutAA
    mut_type = 'R'
    if kinases[acc].fasta[int(uniprot_position)-1] != wtAA:
        raise Exception(f'{kinases[acc].fasta[int(uniprot_position)-1]} found rather than {wtAA} pos {uniprot_position} in {acc}')
    if uniprot_position not in kinases[acc].mutations[mut_type]:
        kinases[acc].mutations[mut_type].append(uniprot_position)

'''ACT/DEACT mutations'''
for line in open('../AK_mut_w_sc_feb2023/act_deact_v2.tsv', 'r'):
    if line.split('\t')[0] == 'uniprot_name':
        continue
    acc = line.split('\t')[1]
    mut_type = line.split('\t')[5]
    wtAA = line.split('\t')[2]
    uniprot_position = line.split('\t')[3]
    if kinases[acc].fasta[int(uniprot_position)-1] != wtAA:
        raise Exception(f'{kinases[acc].fasta[int(uniprot_position)-1]} found rather than {wtAA} pos {uniprot_position} in {acc}')
    if uniprot_position not in kinases[acc].mutations[mut_type]:
        kinases[acc].mutations[mut_type].append(uniprot_position)

dic_mutations = {}
for acc in kinases:
    for mut_type in kinases[acc].mutations:
        if len(kinases[acc].mutations[mut_type]) > 0:
            dic_mutations[acc] = kinases[acc].mutations
            break

with open('kinases_mutations.pkl', 'wb') as f:
    pickle.dump(dic_mutations, f)

with open('kinases_mutations.pkl', 'rb') as f:
    dic_mutations = pickle.load(f)

print (dic_mutations)
