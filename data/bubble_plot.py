#!/usr/bin/env python3
# coding: utf-8
'''
Script to do an hmmsearch of kinases against Pkinase.hmm
and draw a bubble plot for each domain position
'''
import os, sys
from turtle import position
import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact
import plotly.express as px

class Gene:
    '''
    class to define gene and their accessions
    '''
    def __init__(self,) -> None:
        pass
class Mutation:
    '''
    class to define mutations
    '''
    def __init__(self, name, category, num_samples):
        self.name = name
        self.category = category
        self.samples = num_samples
class Kinase:
    '''
    class to define kinases
    '''
    def __init__(self, acc=None, gene=None):
        self.gene = gene
        self.acc = acc
        self.sequence = ''
        self.kinase_to_pfam = {}
        self.act = {}
        self.deact = {}
        self.resistance = {}
    @classmethod
    def from_acc_gene(cls):
        '''
        Alternate way to initialize
        '''
        kinase = cls()
        return kinase
    def get_fasta_formatted(self) -> str:
        return '>'+self.acc+'\n'+self.sequence
    def check_position(self, mutation):
        '''
        Check if the given mutation position and
        WT exist also in the sequence, else raise
        an error
        '''
        wt = mutation[0]
        mutation_position = int(mutation[1:-1])
        for position, aa in enumerate(self.sequence, start=1):
            if position != mutation_position:
                continue
            if aa != wt:
                print (f"position {position} in Gene:{self.gene}, Acc:{self.acc} has {aa}, and not {mutation[0]}; mutation given: {mutation}")
                # raise ValueError(f"position {position} has {aa}, and not {mutation} in {self.gene}, {self.acc}")
                return False
            else:
                return True
    def display(self) -> None:
        '''
        Function to display the instance items
        '''
        print ('kinase:', self.name, sep='\t')
        print ('act:', self.act, sep='\t')
        print ('deact:', self.deact, sep='\t')
        print ('res:', self.resistance, sep='\t')
        print ('kinase_to_pfam', self.kinase_to_pfam)
    
def run_hmmsearch():
    '''
    Fetch all kinase sequences and run HMMSEARCH
    '''
    kinase_dic = {}
    gene_to_accs_dic = {}
    ## Run hmmsearch against Kinase FASTA seqs
    # os.system('hmmsearch -o hmmsearchAlignment3.txt ../pfam/Pkinase.hmm ../KA/fastaForAlignment3.fasta')
    # Fetch all FASTA seqeuences
    fasta_path = '../KA/UniProtFasta2/'
    for file in os.listdir(fasta_path):
        if file.endswith('.fasta') is False:
            continue
        with open(fasta_path+file, 'r') as fp:
            lines = fp.readlines()
            for line in lines:
                if line[0] == '>':
                    acc = line.split('|')[1]
                    gene = line.split('GN=')[1].split()[0]
                    kinase_dic[acc] = Kinase(acc, gene)
                    if gene not in gene_to_accs_dic:
                        gene_to_accs_dic[gene] = [acc]
                    else:
                        if acc not in gene_to_accs_dic[gene]:
                            gene_to_accs_dic[gene].append(acc)
                            gene_to_accs_dic[gene].sort()
                else:
                    kinase_dic[acc].sequence += line.rstrip()
    # Save all sequences together in a FASTA file
    all_kinases_fasta = ''
    for acc in kinase_dic:
        all_kinases_fasta += kinase_dic[acc].get_fasta_formatted() + '\n'
    open('allKinases.fasta', 'w').write(all_kinases_fasta)
    # Run HMMSEARCH against saved sequences
    os.system('hmmsearch -o allKinasesHmmsearch.txt ../pfam/Pkinase.hmm allKinases.fasta')
    return kinase_dic, gene_to_accs_dic

kinase_dic, gene_to_accs_dic = run_hmmsearch()
print (kinase_dic)

## Dictionary that maps gene names to accessions
# acc_to_gene = {}
# for line in open('../../DB/uniprot/uniprot_sprot_human.fasta', 'r'):
#     if line[0] == '>':
#         # print (line)
#         acc = line.split('|')[1]
#         if 'GN=' in line:
#             gene = line.split('GN=')[1].split()[0]
#             acc_to_gene[acc] = gene

## read the HMMSEARCH output
pfam = {}
flag = 0
for line in open('allKinasesHmmsearch.txt', 'r'):
    if len(line.split()) == 0:
        continue
    if line[:2] == '>>':
        ## lines with kinase start
        kinase = line.split('>>')[1].lstrip().rstrip()
        # Raise an error if the kinase instance not found
        if kinase not in kinase_dic:
            raise ValueError(f'{kinase} not found in the HMMSearch output')
        flag = 1
    elif line.split()[0] == 'Pkinase':
        ## lines with Pkinase domain
        pfam_start, pfam_seq, pfam_end = int(line.split()[1]), line.split()[2], int(line.split()[3])
        count = int(line.split()[1])
        for char in pfam_seq:
            if char not in ['.', '-']:
                pfam[count] = char+str(count)
                count += 1
    elif flag == 1:
        if kinase == line.split()[0]:
            ## lines with kinase
            kin_start, kin_seq, kin_end = int(line.split()[1]), line.split()[2], int(line.split()[3])
            for pfam_char, kin_char in zip(pfam_seq, kin_seq):
                if pfam_char not in ['.', '-'] and kin_char not in ['.', '-']:
                    kinase_dic[kinase].kinase_to_pfam[kin_start] = pfam_start
                    pfam_start += 1
                    kin_start += 1
                elif pfam_char in ['.', '-']:
                    kin_start += 1
                elif kin_char in ['.', '-']:
                    pfam_start += 1
                else:
                    print ('Exception found', kinase)
                    sys.exit()
print (kinase_dic['P21802'].kinase_to_pfam[628])
print (kinase_dic['P21802'].kinase_to_pfam[628])
# sys.exit()
## activating/inactivating mutations from UniProt
for line in open('../KA/act_deact_mut_for_scores_fin.tsv', 'r'):
    if line.split('\t')[0] != 'uniprot_id':
        acc = line.split('\t')[0]
        gene = line.split('\t')[5]
        # Raise an error when not found in HMMSearch o/p
        if acc not in kinase_dic:
            raise ValueError(f'{gene} not found in the HMMSearch output')
        # prepare mutation
        mutation = line.split('\t')[1] + line.split('\t')[2] + line.split('\t')[3]
        num_samples = 1
        # Ignore mutations where you see del (deletions)
        if 'del' in mutation:
            continue
        # Check if the given mutation position exist else throw an error
        if kinase_dic[acc].check_position(mutation) == True:
            # Whether activating/deactivating
            status = line.replace('\n', '').split('\t')[-1]
            if status == 'A':
                kinase_dic[acc].act[mutation] = Mutation(mutation, 'activating', num_samples)
            else:
                kinase_dic[acc].deact[mutation] = Mutation(mutation, 'deactivating', num_samples)

## Resistance mutations from COSMIC
'''
for line in open('../KA/resistance_mutations_w_scores_aligned_fin.tsv', 'r'):
    if line.split('\t')[0] != 'Gene.Name':
        gene = line.split('\t')[0]
        accs = gene_to_accs_dic[gene]
        for acc in accs:
            print (accs)
            # Raise an error when not found in HMMSearch o/p
            if acc not in kinase_dic:
                raise ValueError(f'{acc} not found in the HMMSearch output')
            # prepare mutation
            mutation = line.split('\t')[1]
            # Ignore mutations where you see del (deletions)
            if 'del' in mutation:
                continue
            # Check if the given mutation position exist else throw an error
            if kinase_dic[acc].check_position(mutation) == True:
                kinase_dic[acc].resistance[mutation] = Mutation(mutation, 'resistance')
                break
'''
for line in open('../KA/resistant_mutations_Nov22.tsv.gz', 'rt'):
    if line[0] != '#':
        gene = line.split('\t')[0]
        accs = gene_to_accs_dic[gene]
        for acc in accs:
            print (accs)
            # Raise an error when not found in HMMSearch o/p
            if acc not in kinase_dic:
                raise ValueError(f'{acc} not found in the HMMSearch output')
            # prepare mutation
            mutation = line.split('\t')[2]
            num_samples = int(line.split('\t')[4])
            # Ignore mutations where you see del (deletions)
            if 'del' in mutation or '*' in mutation or '_' in mutation or 'dup' in mutation:
                continue
            # Check if the given mutation position exist else throw an error
            if kinase_dic[acc].check_position(mutation) == True:
                kinase_dic[acc].resistance[mutation] = Mutation(mutation, 'resistance', num_samples)
                break

## read the instances
# print (kinase_dic['P25092'].act)
# print (kinase_dic['P25092'].kinase_to_pfam)
data = []
for kinase in kinase_dic:
    # kinase_dic[kinase].display()
    for mutation in kinase_dic[kinase].act:
        mutationInstance = kinase_dic[kinase].act[mutation]
        kin_pos = int(mutation[1:-1])
        if kin_pos in kinase_dic[kinase].kinase_to_pfam:
            pfam_pos = kinase_dic[kinase].kinase_to_pfam[kin_pos]
            pfam_pos = int(pfam_pos)
            num_samples = mutationInstance.samples
            row = []
            # name = kinase if kinase not in acc_to_gene else acc_to_gene[kinase]
            name = kinase_dic[kinase].gene+'/'+kinase
            row.append(name)
            row.append(pfam_pos)
            row.append(pfam[pfam_pos])
            row.append('activating')
            row.append(mutation)
            ## used to display
            row.append(num_samples)
            ## used to calculate size
            row.append(num_samples)
            data.append(row)
    
    for mutation in kinase_dic[kinase].deact:
        mutationInstance = kinase_dic[kinase].deact[mutation]
        kin_pos = int(mutation[1:-1])
        if kin_pos in kinase_dic[kinase].kinase_to_pfam:
            pfam_pos = kinase_dic[kinase].kinase_to_pfam[kin_pos]
            pfam_pos = int(pfam_pos)
            num_samples = mutationInstance.samples
            row = []
            # name = kinase if kinase not in acc_to_gene else acc_to_gene[kinase]
            name = kinase_dic[kinase].gene+'/'+kinase
            row.append(name)
            row.append(pfam_pos)
            row.append(pfam[pfam_pos])
            row.append('deactivating')
            row.append(mutation)
            ## used to display
            row.append(num_samples)
            ## used to calculate size
            row.append(num_samples)
            data.append(row)
    
    for mutation in kinase_dic[kinase].resistance:
        mutationInstance = kinase_dic[kinase].resistance[mutation]
        kin_pos = int(mutation[1:-1])
        if kin_pos in kinase_dic[kinase].kinase_to_pfam:
            pfam_pos = kinase_dic[kinase].kinase_to_pfam[kin_pos]
            pfam_pos = int(pfam_pos)
            num_samples = mutationInstance.samples
            row = []
            # name = kinase if kinase not in acc_to_gene else acc_to_gene[kinase]
            name = kinase_dic[kinase].gene+'/'+kinase
            row.append(name)
            row.append(pfam_pos)
            row.append(pfam[pfam_pos])
            row.append('resistance')
            row.append(mutation)
            ## used to display
            row.append(num_samples)
            ## used to calculate size
            row.append(np.log2(num_samples))
            data.append(row)

## Add ATP binding sites
atp_sites = {}
for line in open('ATP_binding_sites.tsv', 'r'):
    if line[0] == '#':
        continue
    site = int(line.split('\t')[3])
    pdbs = line.split('\t')[6].replace('\n', '').split(';')
    if site not in atp_sites:
        atp_sites[site] = pdbs
    else:
        atp_sites[site] += pdbs
        atp_sites[site] = list(set(atp_sites[site]))

for site in atp_sites:
    row = []
    row.append('ATPbindingSites')
    row.append(int(site))
    row.append(pfam[site])
    row.append('ATP site')
    row.append('ATP site')
    if site in atp_sites:
        row.append(np.log2(len(atp_sites[site])))
        row.append(np.log2(len(atp_sites[site])))
    else:
        row.append(0)
        row.append(0)
    data.append(row)

## Add PTM sites
phospho_sites = {}
acetyl_sites = {}
for line in open('Kinase_psites2.tsv', 'r'):
    if line[0] == '#':
        continue
    pfam_pos = int(line.split('\t')[3])
    # num_ptmsites = int(line.split('\t')[3].replace('\n', ''))
    type_ptm = line.split('\t')[2].split('-')[1]
    dic = phospho_sites if type_ptm == 'p' else acetyl_sites
    if pfam_pos not in dic:
        dic[pfam_pos] = 1
    else:
        dic[pfam_pos] += 1

for count, dic in enumerate([phospho_sites, acetyl_sites]):
    for pfam_pos in dic:
        row = []
        if count == 0:
            row.append('Psites')
        else:
            row.append('Ksites')
        row.append(int(pfam_pos))
        row.append(pfam[pfam_pos])
        if count == 0:
            row.append('P-site')
            row.append('P-site')
            row.append(np.log2(dic[pfam_pos]))
            row.append(np.log2(dic[pfam_pos]))
        else:
            row.append('K-site')
            row.append('K-site')
            row.append(dic[pfam_pos])
            row.append(dic[pfam_pos])
        data.append(row)

df = pd.DataFrame(data=data, columns=['Kinase', 'Pfam_Position', 'Pfam_Residue', 'Category', 'Mutation', 'Num_Samples', 'Num_Samples_Size'])
df = df.sort_values(by=['Kinase'])
print (list(set(df[df['Category']=='activating'].Pfam_Position)))
allRes = list(set(df[df['Category']=='resistance'].Pfam_Position))
allAct = list(set(df[df['Category']=='activating'].Pfam_Position))
allDeact = list(set(df[df['Category']=='deactivating'].Pfam_Position))
resYactY = list(set(allRes).intersection(allAct))
resYdeactY = list(set(allRes).intersection(allDeact))
actYdeactY = list(set(allAct).intersection(allDeact))
resYactYdeactY = list(set(actYdeactY).intersection(allRes))
resYactN = list(set(allRes) - set(allAct))
resYdeactN = list(set(allRes) - set(allDeact))
resNactY = list(set(allAct) - set(allRes))
resNdeactY = list(set(allDeact) - set(allRes))
resNActN = list(set(df.Pfam_Position) - set(allAct) - set(allRes))
print (resYactY, resYactN)
print (resYactYdeactY, 'resYactYdeactY')
print (actYdeactY, 'actYdeactY')
print (resYactY, 'resYactY')
print (resYactN, 'resYactN')
print (resYdeactY, 'resYdeactY')
print (resYdeactN, 'resYdeactN')
oddsratio, pvalue = fisher_exact([[len(resYactY), len(resYdeactY)], [len(resNactY), len(resNdeactY)]])
print (oddsratio, pvalue)
resYactY = set(allRes).intersection(allAct)
ax = sns.scatterplot(data=df, x="Pfam_Position", y="Kinase", hue="Category")
#ax.tick_params(axis='both', which='minor', labelsize=8)
plt.xlabel('Pfam position')
relevant_pfam = list(set(df.Pfam_Position))
#plt.xticks(range(1, len(pfam), 1), range(1, len(pfam), 1), rotation=90, size=5)
plt.grid()
#plt.show()
fig = px.scatter(df, x="Pfam_Position", y="Kinase", color="Category", size="Num_Samples_Size", symbol="Category", hover_data=["Pfam_Residue", "Mutation", "Num_Samples"],
                color_discrete_map={
                "activating": "green",
                "deactivating": "red",
                "resistance": "blue"})
fig.update_yaxes(ticklabelstep=1)
fig.update_layout(
    # hovermode = 'x',
    xaxis = dict(
        tickmode = 'linear',
        tick0 = 1,
        dtick = 2,
        tickfont = dict(
            family = 'Old Standard TT, serif',
            size = 7,
            color = 'black'
            )
    ),
    yaxis = dict(
        tickmode = 'linear',
        tick0 = 1,
        dtick = 1,
        tickfont = dict(
            size = 8.5,
            )
    )
)
fig.update_xaxes(showspikes=True, spikecolor="green", spikethickness=1, spikesnap="cursor", spikemode="across")
fig.update_yaxes(showspikes=True, spikecolor="orange", spikethickness=1)

fig.show()
print (len(pfam))