#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
from tqdm import tqdm

accs = []
for line in open('humanKinases.fasta', 'r'):
    if line.startswith('>') is False: continue
    accs.append(line.strip().split('|')[1])


homologyDir = '/net/home.isilon/ds-russell/mechismoX/analysis/alignments/data/HUMAN/orthologs_only/'
homologyType = ['all_homs', 'orth', 'excl_para', 'spec_para', 'bpso', 'bpsh']
for acc in accs:
    continue
    for homType in homologyType:
        homFile = homologyDir + acc[:4] + '/' + acc + '_' + homType + '.scores.txt.gz'
        if os.path.isfile(homFile) is False:
            print(acc, homType)
            # sys.exit()
        else:
            os.system('cp '+homFile+' homFiles/')

taxonomyDir = '/net/home.isilon/ds-russell/mechismoX/analysis/alignments/data/HUMAN/orthologs_by_taxonomy/'
taxonomyType = ['eukaryotes', 'vertebrates', 'mammals', 'metazoa', 'arthropods']
for acc in tqdm(accs):
    for taxType in taxonomyType:
        if os.path.isdir('taxFiles/'+taxType) is False:
            os.system('mkdir taxFiles/'+taxType)
        taxFile = taxonomyDir + taxType + '/' + acc[:4] + '/' + acc + '.txt'
        if os.path.isfile(taxFile) is False:
            # print(acc, taxType, taxFile)
            taxFile += '.gz'
            # sys.exit()
        os.system('cp ' + taxFile+' taxFiles/'+taxType)

