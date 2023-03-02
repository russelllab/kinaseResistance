#!/usr/bin/env python3
# coding: utf-8

'''
Script to convert CDS 2 protein
'''

import os, sys, gzip

def translate(seq):

    table = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',				
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    protein =""
    if len(seq)%3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            try:
                if codon in table:
                    protein += table[codon]
                else:
                    protein += 'X'
            except:
                print (codon, protein, seq)
                sys.exit()
    else:
        print ('not a multiple of 3', len(seq))
    return protein.replace('_', '')

dic = {}
for line in gzip.open('All_COSMIC_Genes.fasta.gz', 'rt'):
    # print (line.split())
    if line.split() == []:
        continue
    elif line[0] == '>':
        gene_name = line.split('>')[1].split()[0]
        transcript_id = line.split()[1]
        dic[gene_name] = {'transcript': transcript_id, 'sequence': ''}
    else:
        dic[gene_name]['sequence'] += line.replace('\n', '').upper()

out = ''
for line in open('../AK_mut_w_sc_feb2023/res_mut_old.tsv', 'r'):
    if line.split()[0] != 'Gene.Name':
        gene_name = line.split()[0]
        protein_seq = translate(dic[gene_name]['sequence'])
        # print (gene_name, protein_seq)
        out += '>' + gene_name + ':' + dic[gene_name]['transcript'] + '\n'
        out += protein_seq + '\n'
gzip.open('cosmic_kinases.fasta.gz', 'wt').write(out)