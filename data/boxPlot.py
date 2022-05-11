import os, sys

class kinases:
    def __init__(self, name, aln):
        self.name = name
        self.aln = aln
        self.fasta = ''
        self.map = {}
        self.activating = []
        self.deactivating = []

    class mutations:
        def __init__(self, mutation, label):
            self.mutation = mutation
            self.label = [label]

dic = {}
dic['EGFR'] = kinases('EGFR', 'AAABBBCCC')
dic['FGFR'] = kinases('FGFR', 'AAABBBDDD')
inner1 = dic['EGFR'].mutations('A123D', 'resistance')
inner2 = dic['FGFR'].mutations('A456D', 'activating')
print (dic['FGFR'].mutations.mutation)
print (inner2.mutation)

sys.exit()
kin = {}; mut = {}
## Read alignment and map positions
for line in open('../KA/hmmAlignment.aln', 'r'):
    if line.split()!=[] and line[0]!='#' and line.split()!=['//']:
        name = line.split()[0]
        aln = line.split()[1].replace('\n', '')
        if name not in kin:
            kin[name] = kinases(name, aln)
        else:
            kin[name].aln += aln

for name in kin:
    numAln = 0
    numFasta = 0
    for char in kin[name].aln:
        if char != '.' and char != '-':
            kin[name].fasta += char.upper()
            kin[name].map[numFasta] = numAln
            numAln += 1
            numFasta += 1
        else:
            numAln += 1

alnPos = kin['ABL1_ENST00000318560'].map[254]
print (kin['ABL1_ENST00000318560'].map[254])
print (kin['ABL1_ENST00000318560'].aln[alnPos])
print (kin['ABL1'].aln[alnPos])

activatingMutationsAln = []
for line in open('../KA/kinase_activating_mutations_uniprot_with_gene_name.csv', 'r'):
    if line.split(',')[0] != 'Gene':
        #print (line)
        name = line.split(',')[0]
        try:
            fastaPos = int(line.split(',')[3]) - 1
            alnPos = kin[name].map[fastaPos]
            mutation = line.split(',')[2]+line.split(',')[3]+line.split(',')[4]
            activatingMutationsAln.append(alnPos)
            muts = kin[name].mutations(mutation, 'activating')
            kin[name].activating.append(fastaPos)
        except:
            print (line)

activatingMutationsAln = list(set(activatingMutationsAln))

deactivatingMutationsAln = []
for line in open('../KA/kinase_deactivating_mutations_uniprot_with_gene_name.csv', 'r'):
    if line.split(',')[0] != 'Gene':
        #print (line)
        name = line.split(',')[0]
        if name in kin:
            try:
                fastaPos = int(line.split(',')[3]) - 1
                #print (fastaPos)
                alnPos = kin[name].map[fastaPos]
                mutation = line.split(',')[2]+line.split(',')[3]+line.split(',')[4]
                deactivatingMutationsAln.append(alnPos)
                muts = kin[name].mutations(mutation, 'deactivating')
                kin[name].deactivating.append(fastaPos)
            except:
                print (line)

deactivatingMutationsAln = list(set(deactivatingMutationsAln))

intersectionSet = set.intersection(set(activatingMutationsAln), set(deactivatingMutationsAln))

for name in kin:
    for fastaPos in kin[name].deactivating:
        alnPos = kin[name].map[fastaPos]
        if alnPos in deactivatingMutationsAln:
            print (name, kin[name].fasta[fastaPos], fastaPos+1, alnPos+1)

print (kin['EGFR'].map[796])
