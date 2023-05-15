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
import argparse

parser = argparse.ArgumentParser(
                    prog='TrimAli',
                    description='Trims the alignment based on some parameters',
                    epilog='gurdeep')

parser.add_argument('i', help='Input alignment in CLUSTAL format')
parser.add_argument('f', help='path to kinases sequences in FASTA format')
parser.add_argument('s', help='Start of alignment position (32178)')
parser.add_argument('e', help='End of alignment position (32960)')
parser.add_argument('w', help='+/- Window size')
args = parser.parse_args()
INPUT_FILE = args.i
PATH_TO_FASTA = args.f
START_ALN = int(args.s)
END_ALN = int(args.e)
WINDOW = int(args.w)

class Kinase:
    '''A class to store information about a kinase'''
    def __init__(self, acc, real_acc, name, sequence):
        self.acc = acc
        self.real_acc = real_acc
        self.gene = ''
        self.name = name
        self.sequence = sequence
        self.fasta = ''
        self.trimmed_aln_name = None
        self.aln2seq = {}
        self.seq2aln = {}
        self.ptms = {}
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
    region = line.split()[0].split('|')[3]
    acc = line.split()[0].split('|')[1] + '_' + region
    real_acc = line.split()[0].split('|')[1]
    # add sequence to the object
    if acc not in kinases:
        kinases[acc] = Kinase(acc, real_acc, name, sequence)
    else:
        kinases[acc].sequence += sequence

for acc_region in kinases:
    for line in open(PATH_TO_FASTA, 'r'):
        if line.startswith('>'):
            acc = line.split('|')[1]
            gene = line.split('GN=')[1].split()[0]
        
            if acc_region.split('_')[0] == acc:
                kinases[acc_region].gene = gene
                break

for acc_region in kinases:
    if kinases[acc_region].gene == '':
        print ('Gene not found for', acc_region)

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

trimmed_aln_clustal = '# CLUSTAL\n\n'
trimmed_aln_fasta = ''
trimmed_aln_web = '# CLUSTAL\n\n'
for acc, row in zip(aln_accs, df.to_numpy()):
    # print (name, ''.join(row))
    trimmed_aln_clustal += kinases[acc].name + '|'
    trimmed_aln_clustal += kinases[acc].find_fasta_position(START_ALN - WINDOW) + ' '
    kinases[acc].trimmed_aln_name = kinases[acc].name + '|' + kinases[acc].find_fasta_position(START_ALN - WINDOW)
    trimmed_aln_clustal += ''.join(row) + '\n'

    trimmed_aln_fasta += '>' + kinases[acc].name + '|'
    trimmed_aln_fasta += kinases[acc].find_fasta_position(START_ALN - WINDOW) + '\n'
    trimmed_aln_fasta += ''.join(row) + '\n'

    start, end = kinases[acc].name.split('|')[3].split('-')
    num = kinases[acc].find_fasta_position(START_ALN - WINDOW)
    if start == 'start': start = num
    else: start = str(int(start) + int(num) - 1)
    seq = ''.join(row).replace('-', '').replace('.', '')
    end = str(int(start) + len(seq) - 1)
    trimmed_aln_web += kinases[acc].gene + '|' +kinases[acc].real_acc + '|'
    trimmed_aln_web += str(start) + '-' + str(end) + ' '
    trimmed_aln_web += ''.join(row) + '\n'

open(INPUT_FILE.split('.')[0]+'Trimmed.aln', 'w').write(trimmed_aln_clustal)
open(INPUT_FILE.split('.')[0]+'Trimmed.fasta', 'w').write(trimmed_aln_fasta)
open(INPUT_FILE.split('.')[0]+'TrimmedWeb.aln', 'w').write(trimmed_aln_web)

jalview_annotations = 'JALVIEW_ANNOTATION\n'

'''PTM sites'''
PTM_TYPES = []
## Read phospho, acetyl and methyl sites files
for count, files in enumerate(['Phosphorylation_site_dataset.gz',
                                'Acetylation_site_dataset.gz',
                                'Methylation_site_dataset.gz',
                                'Ubiquitination_site_dataset.gz',
                                'Sumoylation_site_dataset.gz',
                                'O-GalNAc_site_dataset.gz',
                                'O-GlcNAc_site_dataset.gz']):
    FLAG = 0
    for line in gzip.open('../data/PSP/'+files, 'rt'):
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
        if '-' in acc:
            continue
        # if acc not in kinases:
        #     continue
        flag = 0
        for accession in kinases:
            if acc == accession.split('_')[0]:
                flag = 1
                break
        if flag == 0: continue

        ptm_site = line.split('\t')[4]
        ptm_position = ptm_site.split('-')[0]
        aa = ptm_position[0]
        ptm_position = ptm_position[1:]
        ptm_type = ptm_site.split('-')[1]
        species = line.split('\t')[6]
        low_throughput = line.split('\t')[10]
        high_throughput = line.split('\t')[11]
        # ignore line when LT_LIT is zero regardless of ptm type
        if low_throughput == '' and high_throughput == '':
            continue
        # ignore line when LT_LIT is < 2 in psites file
        condition_LTLIT = False
        if low_throughput != '':
            condition_LTLIT = False if int(low_throughput) < 2 else True
        condition_MSLIT = False
        if high_throughput != '':
            condition_MSLIT = False if int(high_throughput) < 2 else True
        if condition_LTLIT is not True and condition_MSLIT is not True:
            continue
        # Ignore non-human PTM sites
        if species != 'human':
            continue
        # print (line)
        # print ('--', ptm_position, acc)
        for accession in kinases:
            if acc != accession.split('_')[0]:
                continue
            start = accession.split('_')[1].split('-')[0]
            end = accession.split('_')[1].split('-')[1]
            if start != 'start' and end != 'end':
                if int(ptm_position) >= int(start) and int(ptm_position) <= int(end):
                    ptm_position = int(ptm_position) - int(start) + 1
                    acc = accession
                    break
            elif start != 'start':
                if int(ptm_position) >= int(start):
                    ptm_position = int(ptm_position) - int(start) + 1
                    acc = accession
                    break
            elif end != 'end':
                if int(ptm_position) <= int(end):
                    acc = accession
                    break
            elif start == 'start' and end == 'end':
                acc = accession
                break
        # print (acc, accession, ptm_position)
        if '_' not in acc:
            print (f'PTM position {ptm_position} outside the domain(s) of {acc}')
            continue
        if kinases[acc].fasta[int(ptm_position)-1] != aa:
            # raise Exception(f'{kinases[acc].fasta[int(ptm_position)-1]} found rather than {aa} pos {ptm_position} in {acc}')
            # print (ptm_position, acc)
            print (f'{kinases[acc].fasta[int(ptm_position)-1]} found rather than {aa} pos {ptm_position} in {acc}')
            continue
        if ptm_type not in kinases[acc].ptms:
            kinases[acc].ptms[ptm_type] = []
        if ptm_type not in PTM_TYPES:
            PTM_TYPES.append(ptm_type)
        if ptm_position not in kinases[acc].ptms[ptm_type]:
            kinases[acc].ptms[ptm_type].append(ptm_position)
        aln_position = kinases[acc].seq2aln[int(ptm_position)-1]
        new_aln_position = aln_position-START_ALN+WINDOW + 1
        jalview_annotations += '\t'.join(
                                [
                                'SEQUENCE_GROUP',
                                ptm_type+'-site',
                                str(new_aln_position),
                                str(new_aln_position),
                                '-1',
                                kinases[acc].trimmed_aln_name
                                ]
                                ) + '\n'

# exit()
'''
for line in open('../data/Kinase_psites4.tsv', 'r'):
    if line[0] == '#':
        continue
    acc = line.split('\t')[0]
    if '-' in acc:
        continue
    ptm_site = line.split('\t')[3]
    ptm_position = ptm_site.split('-')[0]
    aa = ptm_position[0]
    ptm_position = ptm_position[1:]
    ptm_type = ptm_site.split('-')[1]
    if kinases[acc].fasta[int(ptm_position)-1] != aa:
        raise Exception(f'{kinases[acc].fasta[int(ptm_position)-1]} found rather than {aa} pos {ptm_position} in {acc}')
    if ptm_type not in kinases[acc].ptms:
        kinases[acc].ptms[ptm_type] = []
    if ptm_type not in PTM_TYPES:
        PTM_TYPES.append(ptm_type)
    if ptm_position not in kinases[acc].ptms[ptm_type]:
        kinases[acc].ptms[ptm_type].append(ptm_position)
    aln_position = kinases[acc].seq2aln[int(ptm_position)-1]
    new_aln_position = aln_position-START_ALN+WINDOW + 1
    jalview_annotations += '\t'.join(
                            [
                            'SEQUENCE_GROUP',
                             ptm_type,
                             str(new_aln_position),
                             str(new_aln_position),
                             '-1',
                             kinases[acc].trimmed_aln_name
                             ]
                            ) + '\n'
'''
# for ptm_type, color in zip(PTM_TYPE, ['', 'Red']):
for ptm_type in PTM_TYPES:
    jalview_annotations += '\t'.join(
                            [
                            'PROPERTIES',
                             ptm_type+'-site',
                             'colour=Yellow',
                             'outlineColour=000000',
                             'displayBoxes=true',
                             'displayText=true',
                             'colourText=false',
                             'showUnconserved=false'
                             ]
                            ) + '\n'

'''Resistance Mutations'''
'''
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
    aln_position = kinases[acc].seq2aln[int(uniprot_position)-1]
    new_aln_position = aln_position-START_ALN+WINDOW + 1
    jalview_annotations += '\t'.join(
                            [
                            'SEQUENCE_GROUP',
                             mut_type,
                             str(new_aln_position),
                             str(new_aln_position),
                             '-1',
                             kinases[acc].trimmed_aln_name
                             ]
                            ) + '\n'
    # break
    # print (new_aln_position, uniprot_mutation, acc, kinases[acc].name)
'''
for line in open('../AK/AK_mut_w_sc_feb2023/res_mut_v3_only_subs_KD_neighb.tsv', 'r'):
    if line.split()[0] == 'uniprot_id':
        continue
    acc = line.split('\t')[0]
    # cosmic_mutation = line.split('\t')[1]
    wtAA = line.split('\t')[1]
    mutAA = line.split('\t')[3]
    uniprot_position = str(line.split('\t')[2].strip())
    uniprot_mutation = wtAA + uniprot_position + mutAA
    mut_type = 'R'
    for accession in kinases:
        if acc != accession.split('_')[0]:
            continue
        start = accession.split('_')[1].split('-')[0]
        end = accession.split('_')[1].split('-')[1]
        if start != 'start' and end != 'end':
            if int(uniprot_position) >= int(start) or int(uniprot_position) <= int(end):
                uniprot_position = int(uniprot_position) - int(start) + 1
                acc = accession
                break
        elif start != 'start':
            if int(uniprot_position) >= int(start):
                uniprot_position = int(uniprot_position) - int(start) + 1
                acc = accession
                break
        elif end != 'end':
            if int(uniprot_position) <= int(end):
                acc = accession
                break
        elif start == 'start' and end == 'end':
            acc = accession
            break
    # print (acc, accession, uniprot_position)
    if '_' not in acc:
        print (f'position {uniprot_position} of {mut_type} {uniprot_mutation} outside the domains of {acc}')
        continue
    if kinases[acc].fasta[int(uniprot_position)-1] != wtAA:
        raise Exception(f'{kinases[acc].fasta[int(uniprot_position)-1]} found rather than {wtAA} pos {uniprot_position} in {acc}')
    if uniprot_position not in kinases[acc].mutations[mut_type]:
        kinases[acc].mutations[mut_type].append(uniprot_position)
    aln_position = kinases[acc].seq2aln[int(uniprot_position)-1]
    new_aln_position = aln_position-START_ALN+WINDOW + 1
    jalview_annotations += '\t'.join(
                            [
                            'SEQUENCE_GROUP',
                             mut_type,
                             str(new_aln_position),
                             str(new_aln_position),
                             '-1',
                             kinases[acc].trimmed_aln_name
                             ]
                            ) + '\n'
    # break
    # print (new_aln_position, uniprot_mutation, acc, kinases[acc].name)

jalview_annotations += '\t'.join(
                            [
                            'PROPERTIES',
                             'R',
                             'colour=Cyan',
                             'outlineColour=000000',
                             'displayBoxes=true',
                             'displayText=true',
                             'colourText=false',
                             'showUnconserved=false'
                             ]
                            ) + '\n'
# open('jalview_annotations2.txt', 'w').write(jalview_annotations)
# exit()

'''ACT/DEACT mutations'''
for line in open('../AK/AK_mut_w_sc_feb2023/act_deact_v2.tsv', 'r'):
    if line.split('\t')[0] == 'uniprot_name':
        continue
    acc = line.split('\t')[1]
    mut_type = line.split('\t')[5]
    wtAA = line.split('\t')[2]
    uniprot_position = line.split('\t')[3]
    uniprot_mutation = line.split('\t')[2]+line.split('\t')[3]+line.split('\t')[4]
    for accession in kinases:
        if acc != accession.split('_')[0]:
            continue
        start = accession.split('_')[1].split('-')[0]
        end = accession.split('_')[1].split('-')[1]
        if start != 'start' and end != 'end':
            if int(uniprot_position) >= int(start) or int(uniprot_position) <= int(end):
                uniprot_position = int(uniprot_position) - int(start) + 1
                acc = accession
                break
        elif start != 'start':
            if int(uniprot_position) >= int(start):
                uniprot_position = int(uniprot_position) - int(start) + 1
                acc = accession
                break
        elif end != 'end':
            if int(uniprot_position) <= int(end):
                acc = accession
                break
        elif start == 'start' and end == 'end':
            acc = accession
            break
    # print (acc, accession, uniprot_position)
    if '_' not in acc:
        print (f'position {uniprot_position} of {mut_type} {uniprot_mutation} outside the domains of {acc}')
        continue
    if kinases[acc].fasta[int(uniprot_position)-1] != wtAA:
        raise Exception(f'{kinases[acc].fasta[int(uniprot_position)-1]} found rather than {wtAA} pos {uniprot_position} in {acc}')
    if uniprot_position not in kinases[acc].mutations[mut_type]:
        kinases[acc].mutations[mut_type].append(uniprot_position)
    aln_position = kinases[acc].seq2aln[int(uniprot_position)-1]
    new_aln_position = aln_position-START_ALN+WINDOW + 1
    jalview_annotations += '\t'.join(
                            [
                            'SEQUENCE_GROUP',
                             mut_type,
                             str(new_aln_position),
                             str(new_aln_position),
                             '-1',
                             kinases[acc].trimmed_aln_name
                             ]
                            ) + '\n'

for mut_type, color in zip(['A', 'D'], ['Green', 'Red']):
    jalview_annotations += '\t'.join(
                            [
                            'PROPERTIES',
                             mut_type,
                             'colour='+color,
                             'outlineColour=000000',
                             'displayBoxes=true',
                             'displayText=true',
                             'colourText=false',
                             'showUnconserved=false'
                             ]
                            ) + '\n'

open('jalview_annotations.txt', 'w').write(jalview_annotations)
exit()
