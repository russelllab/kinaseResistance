#!/usr/bin/env python3
# coding: utf-8
'''
A script to map and plot PTMs on the kinase domain
'''

import sys, gzip
import argparse

parser = argparse.ArgumentParser(description='', epilog='End of help.')
parser.add_argument('d', help='name of domain to search for')
parser.add_argument('h', help='path to HMMSEARCH output')
parser.add_argument('f', help='path to kinases sequences in FASTA format')
parser.add_argument('o', help='path to output file')
args = parser.parse_args()

# set input file to default if not provided
PFAM_DOM = args.d #['humanKinasesHitsSplitTrimmed']
HMMSEARCH_MAPPINGS = args.h # 'humanKinasesHitsSplitHmmsearchTrimmedMappings.tsv.gz'
PATH_TO_FASTA = args.f # 'humanKinases.fasta'
OUT_FILE = args.o #'Kinase_psites_hits_split_trimmed.tsv'
OUT_TEXT = '#Acc\tGene\tPfam_Dom\tPfam_Position\tPfam_Residue\tPTM_type\n'

class Kinases:
    '''
    class to store kinase information
    '''
    def __init__(self, acc, gene) -> None:
        self.acc = acc
        self.gene = gene
        self.sequence = ''
        self.ptms = {}
        self.kinase_to_pfam = {}
    def get_fasta_formatted(self, gene):
        '''
        Function to return FASTA formatted sequence
        '''
        return '>'+self.acc+'\n'+self.sequence

# Fetch sequence of all human kinases
# and save in dictionary with key as gene name
# human_kinases -> gene -> sequence
kinases_dic = {}
human_kinases = {}
for line in open(PATH_TO_FASTA, 'r'):
    if line[0] == '>':
        gene = line.split('GN=')[1].split()[0]
        acc = line.split('|')[1]
        '''
        Just take the gene names, and that is also
        enough to catch the isoform information
        from PSP
        '''
        kinases_dic[acc] = Kinases(acc, gene)
        human_kinases[gene] = ''
    else:
        kinases_dic[acc].sequence += line.strip()
        human_kinases[gene] += line.strip()

# Fetch all PTM sites in kinases and save in dictionary
# with key as accession and value as the class instance
# kinases_dic -> accession -> ptm_type -> [ptm_sites]
for count, row in enumerate([
                                ['Phosphorylation_site_dataset.gz', 'p'],
                                ['Acetylation_site_dataset.gz', 'ac'],
                                ['Methylation_site_dataset.gz', 'me'],
                                ['Ubiquitination_site_dataset.gz', 'ub'],
                                ['Sumoylation_site_dataset.gz', 'sm'],
                                ['O-GalNAc_site_dataset.gz', 'ga'],
                                ['O-GlcNAc_site_dataset.gz', 'gl'],
                                ]):
    FLAG = 0
    files, code = row[0], row[1]
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
        # Ignore non-kinases
        if gene not in human_kinases:
            continue
        species = line.split('\t')[6]
        # Ignore non-human PTM sites
        if species != 'human':
            continue
        acc = line.split('\t')[2]
        if acc not in kinases_dic:
            # print (line.split('\t')[1])
            if 'iso' not in line.split('\t')[1]:
                print (acc)
            continue
        # if acc not in kinases_dic:
        #     accession_address = Kinases(acc, gene, human_kinases[gene])
        #     kinases_dic[acc] = accession_address
        # else:
        accession_address = kinases_dic[acc]
        ptmsite = line.split('\t')[4]
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
        
        if code not in accession_address.ptms:
            accession_address.ptms[code] = []
        accession_address.ptms[code].append(ptmsite)

# Read the HMMSEARCH output and map the PTM sites
# store the results in a dictionary of the class instance
# kinases_dic -> acc -> kinase_to_pfam -> pfam_dom -> uniprot_pos -> domain_pos
pfam = {}
flag = 0
for line in gzip.open(HMMSEARCH_MAPPINGS, 'rt'):
    # print (line)
    if line[0] == '#': continue
    # print (line.split()[0])
    acc = line.split()[0].split('|')[1]
    domainPosition = int(line.split()[4].rstrip())
    domainAA = line.split()[3].rstrip()
    kinasePosition = int(line.split()[2].rstrip())
    kinases_dic[acc].kinase_to_pfam[kinasePosition] = domainPosition
    pfam[domainPosition] = domainAA

    '''if len(line.split()) == 0:
        continue
    if line[:2] == '>>':
        ## lines with kinase start
        kinase = line.split('>>')[1].lstrip().rstrip().split()[0]
        acc = kinase.split('|')[1]
        # Raise an error if the kinase instance not found
        if acc not in kinases_dic:
            # raise ValueError(f'{kinase} not in PSP')
            # print (f'{kinase} not in PSP')
            flag = 0
            continue
        if PFAM_DOM not in kinases_dic[acc].kinase_to_pfam:
            kinases_dic[acc].kinase_to_pfam[PFAM_DOM] = {}
        flag = 1
    elif line.split()[0] == PFAM_DOM and flag == 1:
        ## lines with Pkinase domain
        # print (line)
        pfam_start, pfam_seq, pfam_end = int(line.split()[1]), line.split()[2], int(line.split()[3])
        count = int(line.split()[1])
        for char in pfam_seq:
            if char not in ['.', '-']:
                pfam[count] = char+str(count)
                count += 1
    elif flag == 1:
        if kinase == line.split()[0]:
            # print (acc, kinase, line.split())
            ## lines with kinase
            kin_start, kin_seq, kin_end = line.split()[1], line.split()[2], line.split()[3]
            if kin_start == '-':
                print (line)
                continue
            kin_start, kin_end = int(kin_start), int(kin_end)
            for pfam_char, kin_char in zip(pfam_seq, kin_seq):
                if pfam_char not in ['.', '-'] and kin_char not in ['.', '-']:
                    kinases_dic[acc].kinase_to_pfam[PFAM_DOM][kin_start] = pfam_start
                    pfam_start += 1
                    kin_start += 1
                elif pfam_char in ['.', '-']:
                    kin_start += 1
                elif kin_char in ['.', '-']:
                    pfam_start += 1
                else:
                    print ('Exception found', kinase)
                    sys.exit()'''

# Map PTM sites to the kinase domain and
# For those that don't map, print the domain
# as '-'. It is important to print all PTM sites
# even if they don't map to the domain.
# print (kinases_dic['Q2M2I8'].ptms)
# print (kinases_dic['Q2M2I8'].kinase_to_pfam)
# print (pfam)
for acc in kinases_dic:
    for ptm in kinases_dic[acc].ptms:
        for ptmsite in kinases_dic[acc].ptms[ptm]:
            if kinases_dic[acc].sequence == '':
                print (acc, 'has no sequence')
                continue
            ptmsite_pos = int(ptmsite[1:].split('-')[0])
            if kinases_dic[acc].sequence[ptmsite_pos-1] not in ['S', 'T', 'Y', 'K', 'R']:
                print (acc, 'has', kinases_dic[acc].sequence[ptmsite_pos-1], 'at', ptmsite_pos, 'and not', ptmsite)
                continue
            if ptmsite_pos in kinases_dic[acc].kinase_to_pfam:
                pfam_pos = kinases_dic[acc].kinase_to_pfam[ptmsite_pos]
                # print (acc, pfam_pos)
                OUT_TEXT += acc + '\t' + kinases_dic[acc].gene + '\t' + PFAM_DOM + '\t' + ptmsite + '\t' + str(pfam_pos) + '\t' + str(pfam[pfam_pos]) + '\n'
            else:
                OUT_TEXT += acc + '\t' + kinases_dic[acc].gene + '\t' + PFAM_DOM + '\t' + ptmsite + '\t' + '-' + '\t' + '-' + '\n'

open(OUT_FILE, 'w').write(OUT_TEXT)
print ('Done')