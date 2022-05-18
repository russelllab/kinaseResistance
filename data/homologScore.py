import os, sys

class mutations:
    def __init__(self, acc):
        self.acc = acc
        self.mut = {}

activating = {}
## Activating mutations from UniProt
for line in open('../KA/kinase_activating_mutations_uniprot_with_gene_name.csv', 'r'):
    if line.split(',')[0] != 'Gene':
        acc = line.split(',')[1]
        if acc not in activating:
            activating[acc] = mutations(acc)
        mut = line.split(',')[2]+line.split(',')[3]+line.split(',')[4]
        if mut not in activating[acc].mut:
            activating[acc].mut[mut] = {}
