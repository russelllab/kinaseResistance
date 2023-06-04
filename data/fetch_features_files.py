#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys

accs = []
for line in open('humanKinases.fasta', 'r'):
    if line.startswith('>') is False: continue
    accs.append(line.strip().split('|')[1])


homologyDir = '/net/home.isilon/ds-russell/mechismoX/analysis/alignments/data/HUMAN/orthologs_only/'
homologyType = ['all_homs', 'orth', 'excl_para', 'spec_para', 'bpso', 'bpsh']
structureDir = '/net/home.isilon/ds-russell/mechismoX/analysis/features/data/VLatest/'
structureType = ['iupred', 'dssp-scores', 'mech_intra']
for acc in accs:
    '''
    for homType in homologyType:
        homFile = homologyDir + acc[:4] + '/' + acc + '_' + homType + '.scores.txt.gz'
        if os.path.isfile(homFile) is False:
            print(acc, homType)
            # sys.exit()
        #else:
            #os.system('cp '+homFile+' homFiles/')
    '''
    for strucType in structureType:
        strucFile = structureDir + acc[:4] + '/' + 'AF-' + acc + '-F1-model_v1.' + strucType + '.gz'
        if os.path.isfile(strucFile) is False:
            print (acc, strucType)
        else:
            os.system('cp '+strucFile+' strucFiles/')
