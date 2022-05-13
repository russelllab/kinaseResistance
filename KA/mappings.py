#!/usr/bin/env python3
# coding: utf-8

import os, sys
dic = {}
for line in open('hmmAlignment.aln', 'r'):
    if line.split()!=[]:
        if line[0] != '#' and line.split()[0] != '//':
            name = line.split()[0]
            seq = line.split()[1]
            #print (name, seq)
            if name not in dic:
                dic[name] = ''
            dic[name] += seq

#print (dic['FGFR2_P21802_UniProt'][:1080])
#sys.exit()
l = '#Name\tAlnPosition\tFastaPosition\tFastaAA\n'
for name in dic:
    aln = 1; fasta = 1
    #print (name, dic[name])
    for char in dic[name]:
        #print (name + '\t' + str(aln) + '\t' + str(fasta) + '\t' + str(char) + '\t' + str(len(dic[name].replace('.', '').replace('-',''))) + '\n')
        if char!='.' and char!='-':
            l += str(name) + '\t' + str(aln) + '\t' + str(fasta) + '\t' + str(char) + '\n'
            fasta += 1
        aln += 1

open('hmmAlignmentMappings.tsv', 'w').write(l)
