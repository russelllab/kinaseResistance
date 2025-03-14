# coding: utf-8
#!/usr/bin/env python3

import os, sys, gzip
dic = {}

for line in gzip.open('../data/allKinasesHmmAlign.aln.gz', 'rt'):
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

gzip.open('../data/allKinasesHmmAlignMappings.tsv.gz', 'wt').write(l)
'''
for line in open('muscleAlignment2.clustal', 'r'):
	if line.split()!=[] and line.split()[0] != 'CLUSTAL':
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

open('muscleAlignmentMappings2.tsv', 'w').write(l)
'''
