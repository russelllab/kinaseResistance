#!/usr/bin/env python3

import os, sys, gzip
sys.path.append('../ML')
import fetchData

dic_region = fetchData.getSSdic()
# print (dic_region)

dic_drug_mut = {}
# get drug data from COSMIC
for line in gzip.open('CosmicResistanceMutations.tsv.gz', 'rt'):
    if line.startswith('Sample Name'):
        continue
    gene = line.strip().split('\t')[2]
    if '_' in gene:
        gene = gene.split('_')[0]
    drug = line.strip().split('\t')[5]
    mutation = line.strip().split('\t')[9].split('.')[1]
    if 'del' in mutation or 'ins' in mutation or '?' in mutation or 'du' in mutation or '_' in mutation or len(mutation)<=1:
        continue
    # print (line)
    position = int(mutation[1:-1])
    if gene not in dic_drug_mut:
        dic_drug_mut[gene] = {}
    if position not in dic_drug_mut[gene]:
        dic_drug_mut[gene][position] = []
    dic_drug_mut[gene][position].append(drug)
    dic_drug_mut[gene][position] = list(set(dic_drug_mut[gene][position]))

text = ''
for line in open('atp_table.tsv', 'r'):
    if line.startswith('ligand'):
        text += '\t'.join(line.strip().split('\t')[:4])
        text += '\tRegion\t'
        text += '\t'.join(line.strip().split('\t')[4:])
        text += '\tDrugs\n'
        continue
    line = line.strip().split('\t')
    wtPos = int(line[2])
    gene = line[5]
    if gene not in dic_drug_mut:
        drugs = '-'
    elif wtPos not in dic_drug_mut[gene]:
        drugs = '-'
    else:
        drugs = ';'.join(dic_drug_mut[gene][wtPos])
    alignmentPos = line[3]
    if alignmentPos == '-':
        region = '-'
    else:
        alignmentPos = int(alignmentPos)
        if alignmentPos in dic_region:
            region = ';'.join(dic_region[alignmentPos])
        else:
            region = '-'
    text += '\t'.join(line[:4]) + '\t'
    text += region + '\t'
    text += '\t'.join(line[4:])
    text += '\t' + drugs + '\n'

# print (text)
open('atp_table_MS.tsv', 'w').write(text)