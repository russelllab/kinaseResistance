#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys

accs = []
for line in open('humnankinases.fasta', 'r'):
    if line.startswith('>') is False: continue
    accs.append(line.strip().split('|')[1])


homologyDir = '/net/home.isilon/ds-russell/mechismoX/analysis/alignments/data/HUMAN/orthologs_only/'
homologyType = ['all_homs', 'orth', 'excl_para', 'spec_para', 'bpso', 'bpsh']
for acc in accs:
    for homType in homologyType:
        homFile = homologyDir + acc[:4] + '/' + acc + '_' + homType + '.scores.txt.gz'
        if os.path.isfile(homFile) is False:
            print(acc, homType)
            # sys.exit()