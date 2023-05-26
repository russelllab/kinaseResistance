#!/usr/bin/env python3
# coding: utf-8

'''
A scrip to extract protein kinases that are annotated as
PF00069 and PF07714 in the UniProt database. Fetch their
sequences and save them in a fasta file. Then run the
hmmsearch program to find the domains in the sequences.
Select the one that are above a certain threshold and
save them in a file.
'''

import gzip

acc2pfam = {}
acc = None
with gzip.open('uniprot_sprot_human.dat.gz' ,'rt') as f:
    for line in f:
        # print (line)
        if line.startswith('//'):
            acc = None
            continue
        if line.startswith('AC'):
            if acc is None:
                acc = line.split()[1].rstrip(';')
        elif line.startswith('DR') and acc is not None:
            # if acc == 'O00141': print (line)
            if 'Pfam' in line:
                pfam = line.split()[2].rstrip(';')
                if pfam == 'PF00069' or pfam == 'PF07714':
                    acc2pfam[acc] = pfam
                    # print (acc, pfam)
                    continue

# print ((acc2pfam)['O00141'])
for acc in acc2pfam:
    print (acc, sep='|')