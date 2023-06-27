#!/usr/bin/env python3

'''
Script to make secondary structure annotations for an alignment
BRAF_HUMAN is used as the reference sequence
'''

import gzip
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('alignment', help='alignment file')
args = parser.parse_args()

alignment_file = args.alignment
output_file = alignment_file.split('.')[0]
if '/' in output_file:
    output_file = output_file.split('/')[-1]
output_file += '_ss.tsv.gz'

ALLOWED_REGIONS = ["N-term-β1","β1","Gly-rich-loop",
                   "β1/β2-loop","β2","β2/β3-loop",
                   "β3","Catalytic-Lys","β3/αC-loop",
                   "αC","αC/β4-loop","β4",
                   "β4/β5-loop","β5","β5/αD-loop",
                   "αD","αD/αE-loop","αE",
                   "Catalytic-loop","HrD-motif",
                   "A-loop","DFG-motif","APE-motif",
                   "αF","αF/αG-loop","αG","αG/αH-loop",
                   "αH","αH/αI-loop","αI","C-term-αI"]

# read ss definitions of BRAF_HUMAN
dic_ss_braf = {}
for line in open('../data/BRAF_ss.tsv', 'r', encoding='utf-8'):
    if line.startswith('#'):
        continue
    region = line.split()[0]
    if region not in ALLOWED_REGIONS:
        raise ValueError('Unexpected region -' + region + '- found in BRAF_ss.tsv')
    if '-' in line.split()[1]:
        start, end = line.split()[1].split('-')
    else:
        start, end = int(line.split()[1]), int(line.split()[1])
    
    for position in range(int(start), int(end)+1):
        if position not in dic_ss_braf:
            dic_ss_braf[position] = []
        dic_ss_braf[position].append(region)

# read alignment
BRAF_ALN = ''
FLAG = 0
for line in open(alignment_file, 'r', encoding='utf-8'):
    if line.startswith('>'):
        # When done reading BRAF_HUMAN, break
        if FLAG == 1:
            break
        name = line.lstrip('>').rstrip('\n')
        # When BRAF_HUMAN is reached, start reading
        if 'BRAF_HUMAN' in name:
            FLAG += 1
            braf_start = int(name.split('|')[4])
            continue
    if FLAG != 1: continue
    BRAF_ALN += line.rstrip('\n').upper()

# print (braf_start, BRAF_ALN)
COUNTER = 0
dic_ss_aln = {}
for aln_pos, char in enumerate(BRAF_ALN):
    if BRAF_ALN[aln_pos] in ['-', '.']:
        continue
    braf_pos = braf_start + COUNTER
    braf_ss = dic_ss_braf[braf_pos]
    for region in braf_ss:
        if region not in dic_ss_aln:
            dic_ss_aln[region] = []
        dic_ss_aln[region].append(aln_pos+1)
        # print (braf_pos, dic_ss_braf[braf_pos], sep='\t')
    COUNTER += 1

# print (dic_ss_aln)

OUT_TEXT = '#region\tstart-end\n'
for region in dic_ss_aln:
    start_region = min(dic_ss_aln[region])
    end_region = max(dic_ss_aln[region])
    if '-loop' in region and '/' in region:
        nterm_region = region.split('-')[0].split('/')[0]
        cterm_region = region.split('-')[0].split('/')[1]
        end_nterm_region = max(dic_ss_aln[nterm_region])
        start_cterm_region = min(dic_ss_aln[cterm_region])
        start_region = end_nterm_region + 1
        end_region = start_cterm_region - 1
    elif region in ['N-term-β1']:
        cterm_region = region.split('-')[2]
        start_cterm_region = min(dic_ss_aln[cterm_region])
        end_region = start_cterm_region - 1
        start_region = 1
    elif region in ['C-term-αI']:
        nterm_region = region.split('-')[2]
        end_nterm_region = max(dic_ss_aln[nterm_region])
        start_region = end_nterm_region + 1
        end_region = len(BRAF_ALN)
    print (region, str(start_region)+'-'+str(end_region), sep='\t')
    OUT_TEXT += region + '\t' + str(start_region) + '-' + str(end_region) + '\n'

gzip.open(output_file, 'wt', encoding='utf-8').write(OUT_TEXT)