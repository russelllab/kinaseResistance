'''
Script to map UniProt seq positions
to Pkinase domain positions
'''
import os, gzip, sys
import argparse, threading
import os, sys, gzip

flag = 0
dic_kinase = {}
for line in gzip.open('humanKinasesHmmsearch2.txt.gz', 'rt'):
	if len(line.split()) == 0:
		continue
	if 'Domain annotation for each sequence (and alignments):' in line:
		flag = 1
		continue
	if flag == 0:
		continue
	if line[:2] == '>>':
		acc = line.split('>>')[1].lstrip().rstrip().rstrip('\n')
		acc = acc.split()[0]
		if acc not in dic_kinase: dic_kinase[acc] = ''
		# if len(dic_kinase) == 2: break
		# print (acc)
		continue
	if line.split()[0] == 'Pkinase':
		start_pfam = int(line.split()[1])
		seq_pfam = line.split()[2]
		end_pfam = int(line.split()[3].rstrip('\n'))
		continue
	if line.split()[0] == acc:
		# print (line)
		start_acc = int(line.split()[1])
		seq_acc = line.split()[2]
		end_acc = int(line.split()[3].rstrip('\n'))
		for pfam_aa, acc_aa in zip(seq_pfam, seq_acc):
			if pfam_aa not in ['.', '-'] and acc_aa not in ['.', '-']:
				dic_kinase[acc] += acc + '\t'
				dic_kinase[acc] += str(acc_aa) + '\t' + str(start_acc) + '\t'
				dic_kinase[acc] += str(pfam_aa) + '\t' + str(start_pfam) + '\n'
				start_acc += 1
				start_pfam += 1
			elif pfam_aa not in ['.', '-']:
				start_pfam += 1
			elif acc_aa not in ['.', '-']:
				start_acc += 1
		# print (dic_kinase)
		continue
# print (text)
text = '#Acc\tUniProt_AA\tUniProt_Pos\tPfam_AA\tPfam_Pos\n'
for acc in dic_kinase:
	# print (dic_kinase[acc])
	text += dic_kinase[acc]

gzip.open('humanKinasesHmmsearchMappings2.tsv.gz', 'wt').write(text)