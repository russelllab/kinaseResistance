#!/usr/bin/env python3

''' A script to take an alignment as input
and trim based on some paramerters'''

import pandas as pd
import numpy as np
from tqdm import tqdm

INPUT_FILE = 'humanKinasesHITS.clustal'
GAP_THRESHOLD = 0.9

class Kinase:
    '''A class to store information about a kinase'''
    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence
        self.fasta = ''
        self.aln2seq = {}
        self.seq2aln = {}
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
    # add sequence to the object
    if name not in kinases:
        kinases[name] = Kinase(name, sequence)
    else:
        kinases[name].sequence += sequence

# convert the alignment to fasta
# and make the mapping dictionaries
aln_data = []
aln_names = []
for num, name in enumerate(kinases):
    # print (kinases['sp|Q92772|CDKL2_HUMAN'].fasta)
    # print (name)
    kinases[name].convert2fasta()
    kinases[name].make_aln2seq()
    aln_data.append(list(kinases[name].sequence))
    aln_names.append(name)
    # if num == 2: break
    
df = pd.DataFrame(aln_data, index=aln_names)
for col in tqdm(df.columns):
    # print (df[col].to_numpy())
    column = df[col].to_numpy()
    sum_gaps = np.count_nonzero(column == '-') + np.count_nonzero(column == '.')
    total = len(column)
    if float(sum_gaps)/total >= GAP_THRESHOLD:
        df = df.drop(columns=col)

# print (df)
df.to_csv('trimmedCols.tsv', sep='\t')

trim2 = ''
for line in open('trimmedCols.tsv', 'r'):
    if line.split('\t')[0] == '':
        continue
    name = line.split('\t')[0]
    aln_seq =  ''.join(line.split('\t')[1:])
    trim2 += '>' + name + '\n' + aln_seq

open('trim2.fasta', 'w').write(trim2)