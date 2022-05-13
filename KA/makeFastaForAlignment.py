import os, sys

l = ''
for line in open('fasta_cosmic_kinases_for_alignment_full.fasta', 'r'):
    l += line

for files in os.listdir('UniProtFasta'):
    if files.endswith('_UniProt.fasta'):
        for line in open('UniProtFasta/'+files, 'r'):
            l += line
        l += '\n'

open('fastaForAlignment.fasta', 'w').write(l)
