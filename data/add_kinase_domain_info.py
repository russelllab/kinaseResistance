#!/usr/bin/env python3

import os, sys, gzip

def runHmmsearch(acc):
    pdomains = ['Pkinase', 'PK_Tyr_Ser-Thr']
    row = []
    for pdomain in pdomains:
        os.system('hmmsearch --tblout hh '+ ' ../pfam/' + pdomain +'.hmm '\
                  '../KA/UniProtFasta2/'+ acc + '.fasta.gz')
        for line in open('hh', 'r'):
            if line.startswith('#'): continue
            if line.split()[0] == acc: continue
            print (line.split())
            # print (line.split('\t')[4], acc)
            if float(line.split()[4]) > 1e-5: continue
            row.append(line.split()[3])
    if len(row) == 0:
        print (acc, 'nothing works')
        sys.exit()
    else:
        return row

path2fasta = sys.argv[1]
dic_kinase = {}
for files in os.listdir(path2fasta):
    if files.endswith('.txt.gz') == False: continue
    acc = files.split('.')[0]
    if '-' in acc: continue
    if acc not in dic_kinase: dic_kinase[acc] = []
    for line in gzip.open(path2fasta + files, 'rt'):
        if line.startswith('DR') == False: continue
        if line.split(';')[0].split()[1] != 'Pfam': continue
        if line.split(';')[1].split()[0] in ['PF00069', 'PF07714', 'PF00454',\
                                             'PF02518', 'PF02816', 'PF06743',\
                                            'PF01163', 'PF10494', 'PF03109']:
            dic_kinase[acc].append(line.split(';')[1].split()[0])

for acc in dic_kinase:
    if len(dic_kinase[acc]) == 0:
        if acc in ['O60885', 'Q9UIG0', 'P21675', 'Q13263',\
                   'O15164', 'Q9NRL2', 'P53004', 'Q9Y5P4',\
                    'Q5VZY9', 'Q9UPN9', 'P11274', 'Q8NI60',\
                    'Q12979', 'Q12979', 'Q15059', 'Q58F21',\
                    'P25440', 'Q58F21', 'Q8IZX4']:
            dic_kinase[acc] = '-'
        else:
            dic_kinase[acc] = runHmmsearch(acc)