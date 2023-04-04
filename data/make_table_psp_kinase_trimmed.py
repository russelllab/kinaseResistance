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
# PFAM_DOMS = ['humanKinasesTrimmed'] # Domains to search for
PFAM_DOMS = ['humanKinasesHitsSplitTrimmed'] # Domains to search for
OUT_TEXT = '#Acc\tGene\tPfam_Dom\tPfam_Position\tPfam_Residue\tPTM_type\n'
OUT_FILE = 'Kinase_psites_hits_split_trimmed.tsv'

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
        self.msites = []
        self.usites = []
        self.sumosites = []
        self.gasites = []
        self.glsites = []
        self.kinase_to_pfam = {}
    def get_fasta_formatted(self, gene):
        '''
        Function to return FASTA formatted sequence
        '''
        return '>'+self.acc+'\n'+self.sequence

## fetch all human kinases
human_kinases = []
for line in open('humanKinases.fasta', 'r'):
    if line[0] != '>':
        continue
    gene = line.split('GN=')[1].split()[0]
    '''
    Just take the gene names, and that is also
    enough to catch the isoform information
    from PSP
    '''
    human_kinases.append(gene)

## fetch all PTM sites in kinases
kinases_dic = {}
## Read phospho, acetyl and methyl sites files
for count, files in enumerate(['Phosphorylation_site_dataset.gz',
                                'Acetylation_site_dataset.gz',
                                'Methylation_site_dataset.gz',
                                'Ubiquitination_site_dataset.gz',
                                'Sumoylation_site_dataset.gz',
                                'O-GalNAc_site_dataset.gz',
                                'O-GlcNAc_site_dataset.gz']):
    FLAG = 0
    for line in gzip.open('PSP/'+files, 'rt'):
        if len(line.split()) == 0:
            continue
        # flag when header found
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
        high_throughput = line.split('\t')[11]
        # ignore line when LT_LIT is zero regardless of ptm type
        if low_throughput == '' and high_throughput == '':
            continue
        # ignore line when LT_LIT is < 1 in psites file
        condition_LTLIT = False
        if low_throughput != '':
            condition_LTLIT = False if int(low_throughput) < 1 else True
        condition_MSLIT = False
        if high_throughput != '':
            condition_MSLIT = False if int(high_throughput) < 2 else True
        if condition_LTLIT is not True and condition_MSLIT is not True:
            continue
        # Ignore non-human PTM sites
        if species != 'human':
            continue
        # Ignore non-kinases
        if gene not in human_kinases:
            continue

        if acc not in kinases_dic:
            accession_address = Kinases(acc, gene)
            kinases_dic[acc] = accession_address
        else:
            accession_address = kinases_dic[acc]
        
        if count == 0:
            accession_address.psites.append(ptmsite)
        elif count == 1:
            accession_address.ksites.append(ptmsite)
        elif count == 2:
            accession_address.msites.append(ptmsite)
        elif count == 3:
            accession_address.usites.append(ptmsite)
        elif count == 4:
            accession_address.sumosites.append(ptmsite)
        elif count == 5:
            accession_address.gasites.append(ptmsite)
        elif count == 6:
            accession_address.glsites.append(ptmsite)
# print (len(kinases_dic.keys()))
# print (kinases_dic.keys())
# sys.exit()

## retreive and save all kinase
## sequences in FASTA format
all_kinases_fasta = ''
for acc in kinases_dic:
    if acc in ['Q59GL6', 'Q15300', 'Q8IWY7']:
        continue
    if os.path.isfile(PATH_TO_FASTA+acc+'.fasta') is False and os.path.isfile(PATH_TO_FASTA+acc+'.fasta.gz') is False:
        os.system('wget -O '+PATH_TO_FASTA+acc+'.fasta '+'https://rest.uniprot.org/uniprotkb/'+acc+'.fasta')
    with gzip.open(PATH_TO_FASTA+acc+'.fasta.gz', 'rt') as fp:
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
'''
open('allKinases.fasta', 'w').write(all_kinases_fasta)
'''
# Run HMMSEARCH against saved sequences
os.system('hmmsearch -o allKinasesHitsSplitHmmsearchTrimmed.txt ../pfam/humanKinasesHitsSplitTrimmed.hmm allKinases.fasta')

# print(kinases_dic)
## read the HMMSEARCH output
for pfam_dom in PFAM_DOMS:
    HMMSEARCH_OUT = 'allKinasesHitsSplitHmmsearchTrimmed.txt'
    pfam = {}
    flag = 0
    for line in open(HMMSEARCH_OUT, 'r'):
        # print (line)
        if len(line.split()) == 0:
            continue
        if line[:2] == '>>':
            ## lines with kinase start
            kinase = line.split('>>')[1].lstrip().rstrip()
            # Raise an error if the kinase instance not found
            if kinase not in kinases_dic:
                # raise ValueError(f'{kinase} not in PSP')
                # print (f'{kinase} not in PSP')
                flag = 0
                continue
            if pfam_dom not in kinases_dic[kinase].kinase_to_pfam:
                kinases_dic[kinase].kinase_to_pfam[pfam_dom] = {}
            flag = 1
        elif line.split()[0] == pfam_dom and flag == 1:
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
                kin_start, kin_seq, kin_end = line.split()[1], line.split()[2], line.split()[3]
                if kin_start == '-':
                    print (line)
                    continue
                kin_start, kin_end = int(kin_start), int(kin_end)
                for pfam_char, kin_char in zip(pfam_seq, kin_seq):
                    if pfam_char not in ['.', '-'] and kin_char not in ['.', '-']:
                        kinases_dic[kinase].kinase_to_pfam[pfam_dom][kin_start] = pfam_start
                        pfam_start += 1
                        kin_start += 1
                    elif pfam_char in ['.', '-']:
                        kin_start += 1
                    elif kin_char in ['.', '-']:
                        pfam_start += 1
                    else:
                        print ('Exception found', kinase)
                        sys.exit()

    ## Map PTM sites to kinases domain and
    ## save in corresponding dictionaries
    pfam_phospho = {}
    pfam_acetyl = {}
    pfam_methyl = {}
    pfam_ubiq = {}
    pfam_sumo = {}
    pfam_ga = {}
    pfam_gl = {}
    for count, dic in enumerate([pfam_phospho, pfam_acetyl, pfam_methyl, pfam_ubiq, pfam_sumo, pfam_ga, pfam_gl]):
        for acc in kinases_dic:
            if pfam_dom not in kinases_dic[acc].kinase_to_pfam:
                    continue
            if count == 0:
                ptm_sites = kinases_dic[acc].psites
            elif count == 1:
                ptm_sites = kinases_dic[acc].ksites
            elif count == 2:
                ptm_sites = kinases_dic[acc].msites
            elif count == 3:
                ptm_sites = kinases_dic[acc].usites
            elif count == 4:
                ptm_sites = kinases_dic[acc].sumosites
            elif count == 5:
                ptm_sites = kinases_dic[acc].gasites
            elif count == 6:
                ptm_sites = kinases_dic[acc].glsites
            for ptmsite in ptm_sites:
                if kinases_dic[acc].sequence == '':
                    continue
                # print (acc, psite)
                ptmsite_pos = int(ptmsite[1:].split('-')[0])
                # print (ptmsite_pos)
                if kinases_dic[acc].sequence[ptmsite_pos-1] not in ['S', 'T', 'Y', 'K', 'R']:
                    print (acc, 'has', kinases_dic[acc].sequence[ptmsite_pos-1], 'at', ptmsite_pos, 'and not', ptmsite)
                    # sys.exit()
                    continue
                if ptmsite_pos in kinases_dic[acc].kinase_to_pfam[pfam_dom]:
                    # print (acc, psite, psite_pos)
                    pfam_pos = kinases_dic[acc].kinase_to_pfam[pfam_dom][ptmsite_pos]
                    OUT_TEXT += acc + '\t' + kinases_dic[acc].gene + '\t' + pfam_dom + '\t' + ptmsite + '\t' + str(pfam_pos) + '\t' + str(pfam[pfam_pos]) + '\n'
                    if pfam_pos not in dic:
                        dic[pfam_pos] = 1 
                    else:
                        dic[pfam_pos] += 1

open(OUT_FILE, 'w').write(OUT_TEXT)
print ('Done')