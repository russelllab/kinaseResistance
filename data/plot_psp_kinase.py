#!/usr/bin/env python3
# coding: utf-8
'''
A script to map and plot PTMs on the kinase domain
'''

import os, sys, gzip
import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact
import plotly.express as px

PATH_TO_FASTA = '../KA/UniProtFasta2/'

class Kinases:
    '''
    class to define accessions
    '''
    def __init__(self, acc, gene) -> None:
        self.acc = acc
        self.gene = gene
        self.sequence = ''
        self.psites = []
        self.ksites = []
        self.kinase_to_pfam = {}
    def get_fasta_formatted(self, gene):
        '''
        Fucntion to return FASTA formatted sequence
        '''
        return '>'+self.acc+'\n'+self.sequence

## fetch all human kinases
human_kinases = []
for line in open('humanKinases.fasta', 'r'):
    if line[0] != '>':
        continue
    gene = line.split('GN=')[1].split()[0]
    human_kinases.append(gene)

## fetch all PTM sites in kinases
kinases_dic = {}
for count, files in enumerate(['Phosphorylation_site_dataset.gz', 'Acetylation_site_dataset.gz']):
    FLAG = 0
    for line in gzip.open('PSP/'+files, 'rt'):
        if len(line.split()) == 0:
            continue
        if line.split()[0] == 'GENE':
            FLAG = 1
            continue
        if FLAG != 1:
            continue
        gene = line.split('\t')[0]
        # print (line)
        acc = line.split('\t')[2]
        ptmsite = line.split('\t')[4]
        species = line.split('\t')[6]
        low_throughput = line.split('\t')[10]
        if low_throughput == '':
            continue
        if int(low_throughput) < 2 and count == 0:
            continue
        if int(low_throughput) < 0 and count == 1:
            continue
        if species != 'human':
            continue
        if gene not in human_kinases:
            continue
        if acc not in kinases_dic:
            accession_address = Kinases(acc, gene)
            kinases_dic[acc] = accession_address
        else:
            accession_address = kinases_dic[acc]
        
        if count == 0:
            accession_address.psites.append(ptmsite)
        else:
            accession_address.ksites.append(ptmsite)
# print (len(kinases_dic.keys()))
# print (kinases_dic.keys())
# sys.exit()

## get all kinases' FASTA formatted sequence
all_kinases_fasta = ''
for acc in kinases_dic:
    if acc in ['Q59GL6', 'Q15300', 'Q8IWY7']:
        continue
    if os.path.isfile(PATH_TO_FASTA+acc+'.fasta') is False:
        os.system('wget -O '+PATH_TO_FASTA+acc+'.fasta '+'https://rest.uniprot.org/uniprotkb/'+acc+'.fasta')
    with open(PATH_TO_FASTA+acc+'.fasta', 'r') as fp:
        lines = fp.readlines()
        for line in lines:
            if line[0] == '>':
                acc = line.split('|')[1]
                gene = line.split('GN=')[1].split()[0]
                # print (acc, gene)
                accession_address = kinases_dic[acc]
            else:
                accession_address.sequence += line.replace('\n', '')
        all_kinases_fasta += accession_address.get_fasta_formatted(gene) + '\n'

open('allKinases.fasta', 'w').write(all_kinases_fasta)

# Run HMMSEARCH against saved sequences
os.system('hmmsearch -o allKinasesHmmsearch.txt ../pfam/Pkinase.hmm allKinases.fasta')

## read the HMMSEARCH output
pfam = {}
flag = 0
for line in open('allKinasesHmmsearch.txt', 'r'):
    if len(line.split()) == 0:
        continue
    if line[:2] == '>>':
        ## lines with kinase start
        kinase = line.split('>>')[1].lstrip().rstrip()
        # Raise an error if the kinase instance not found
        if kinase not in kinases_dic:
            raise ValueError(f'{kinase} not found in the HMMSearch output')
        flag = 1
    elif line.split()[0] == 'Pkinase':
        ## lines with Pkinase domain
        pfam_start, pfam_seq, pfam_end = int(line.split()[1]), line.split()[2], int(line.split()[3])
        count = int(line.split()[1])
        for char in pfam_seq:
            if char not in ['.', '-']:
                pfam[count] = char+str(count)
                count += 1
    elif flag == 1:
        if kinase == line.split()[0]:
            ## lines with kinase
            kin_start, kin_seq, kin_end = int(line.split()[1]), line.split()[2], int(line.split()[3])
            for pfam_char, kin_char in zip(pfam_seq, kin_seq):
                if pfam_char not in ['.', '-'] and kin_char not in ['.', '-']:
                    kinases_dic[kinase].kinase_to_pfam[kin_start] = pfam_start
                    pfam_start += 1
                    kin_start += 1
                elif pfam_char in ['.', '-']:
                    kin_start += 1
                elif kin_char in ['.', '-']:
                    pfam_start += 1
                else:
                    print ('Exception found', kinase)
                    sys.exit()

pfam_phospho = {}
pfam_acetyl = {}
for count, dic in enumerate([pfam_phospho, pfam_acetyl]):
    for acc in kinases_dic:
        if count == 0:
            ptm_sites = kinases_dic[acc].psites
        else:
            ptm_sites = kinases_dic[acc].ksites
        for ptmsite in ptm_sites:
            if kinases_dic[acc].sequence == '':
                continue
            # print (acc, psite)
            ptmsite_pos = int(ptmsite[1:].split('-')[0])
            if kinases_dic[acc].sequence[ptmsite_pos-1] not in ['S', 'T', 'Y', 'K']:
                print (acc, 'has', kinases_dic[acc].sequence[ptmsite_pos-1], 'at', ptmsite_pos, 'and not', ptmsite)
                # sys.exit()
                continue
            if ptmsite_pos in kinases_dic[acc].kinase_to_pfam:
                # print (acc, psite, psite_pos)
                pfam_pos = kinases_dic[acc].kinase_to_pfam[ptmsite_pos]
                if pfam_pos not in dic:
                    dic[pfam_pos] = 1 
                else:
                    dic[pfam_pos] += 1

# print (pfam_phospho)
# print (pfam_acetyl)
# sys.exit()
out_text = '#Pfam_Position\tPfam_Residue\tType_PTM\tNum_PTMsites\n'
for i in range(1, 265):
    for count, dic in enumerate([pfam_phospho, pfam_acetyl]):
        if i in dic:
            if count == 0:
                type_ptm = 'p-site'
            else:
                type_ptm = 'k-site'
            print (i, pfam[i], type_ptm, dic[i])
            out_text += str(i) + '\t' + pfam[i] + '\t' + type_ptm + '\t' + str(dic[i]) + '\n'

open('Kinase_psites.tsv', 'w').write(out_text)