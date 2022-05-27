import os, sys, requests

accs = []
for line in open('humanKinases.txt', 'r'):
    accs.append(line.replace('\n',''))

fasta = {}; flag = 0
for line in open('/home/gurdeep/projects/DB/uniprot/uniprot_sprot_human.fasta', 'r'):
    if line[0] == '>':
        acc = line.split('|')[1]
        flag = 0
        if acc in accs:
            fasta[acc] = line
            flag = 1
    elif flag == 1:
        fasta[acc] += line

print (len(list(set(accs))))
print (len(fasta))

l = ''
for acc in fasta:
    l += fasta[acc]

open('humanKinases.fasta', 'w').write(l)
