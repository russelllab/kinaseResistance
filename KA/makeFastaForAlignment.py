import os, sys

l = ''
'''
for line in open('fasta_cosmic_kinases_for_alignment_full.fasta', 'r'):
    l += line
'''
for files in os.listdir('UniProtFasta2'):
    #if files.endswith('_UniProt.fasta'):
    if files.endswith('.fasta'):
        for line in open('UniProtFasta2/'+files, 'r'):
            l += line
        l += '\n'

open('fastaForAlignment2.fasta', 'w').write(l)
