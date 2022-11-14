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

PFAM_DOM = 'PK_Tyr_Ser-Thr'
# PFAM_DOM = 'Pkinase'
HMMSEARCH_OUT = 'allKinasesHmmsearch'+PFAM_DOM+'.txt'

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
    def __init__(self, acc=None, gene=None, pfam_domains=None):
        self.gene = gene
        self.acc = acc
        self.pfam_domains = pfam_domains
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
        acc = file.split('.')[0]
        # print (acc+'.txt')
        if os.path.isfile(fasta_path+acc+'.txt') is False:
            os.system('wget -O ' + fasta_path + acc+'.txt ' + 'https://rest.uniprot.org/uniprotkb/'+acc+'.txt')
        pfam_domains = []
        canonical_acc = acc if '-' not in acc else acc.split('-')[0]
        # print (canonical_acc, acc.split(), file, fasta_path+acc+'.txt')
        for line in open(fasta_path + canonical_acc + '.txt', 'r'):
            if line[:2] != 'DR':
                continue
            if 'Pfam;' in line.split('DR')[1].split():
                # print (line.split('DR')[1].split())
                pfam_domains.append(line.split('DR')[1].split()[2].replace(';', ''))
        if len(pfam_domains) == 0:
            for line in open('allKinasesHmmsearch'+PFAM_DOM+'.txt', 'r'):
                if line[:2] == '>>':
                    protein = line.split('>>')[1].lstrip().rstrip().replace('\n', '')
                    if protein in [acc, canonical_acc]:
                        pfam_domains = [PFAM_DOM]
                        break
            print (acc, canonical_acc, 'found no PFAM accession')
            continue
        with open(fasta_path+file, 'r') as fp:
            lines = fp.readlines()
            for line in lines:
                if line[0] == '>':
                    acc = line.split('|')[1]
                    gene = line.split('GN=')[1].split()[0]
                    kinase_dic[acc] = Kinase(acc, gene, pfam_domains)
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
    os.system('hmmsearch -o ' + HMMSEARCH_OUT + ' ../pfam/' + PFAM_DOM + '.hmm' + ' allKinases.fasta')
    return kinase_dic, gene_to_accs_dic

kinase_dic, gene_to_accs_dic = run_hmmsearch()
# print (kinase_dic['P36888'].pfam_domains)
# sys.exit()
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
for line in open(HMMSEARCH_OUT, 'r'):
    if len(line.split()) == 0:
        continue
    if line[:2] == '>>':
        ## lines with kinase start
        kinase = line.split('>>')[1].lstrip().rstrip()
        # Raise an error if the kinase instance not found
        if kinase not in kinase_dic:
            raise ValueError(f'{kinase} not found in the HMMSearch output')
        flag = 1
    elif line.split()[0] == PFAM_DOM:
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
        num_samples = len(line.split('\t')[4].split('PubMed:')) - 1
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
# print (kinase_dic['Q9UM73'].resistance['C1156Y'].samples)
# print (kinase_dic['P25092'].kinase_to_pfam)
data = []
for kinase in kinase_dic:
    if PFAM_DOM not in kinase_dic[kinase].pfam_domains:
        # print (PFAM_DOM, kinase_dic[kinase].pfam_domains)
        continue
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

## ATP binding sites
ligand_sites = {}
LIGANDS = ['ATP', 'ADP', 'MG', 'MN',
            '0WM','1LT','07J','DB8',
            '6GY','4MK','6T2','VGH',
            'P06','1N1','AQ4','E53',
            'IRE','STI','NIL','YY3',
            'LQQ','P30','BAX','B49',
            '032']
for line in open('ATP_binding_sites3.tsv', 'r'):
    if line[0] == '#':
        continue
    # Consider only the one with the specified domain
    if line.split('\t')[2] != PFAM_DOM:
        continue
    ligand = line.split('\t')[3]
    # Mark ADP as ATP and MN as MG
    if ligand == 'ADP':
        ligand = 'ATP'
    elif ligand == 'MN':
        ligand = 'MG'
    elif ligand not in ['ATP', 'ADP', 'MG', 'MN']:
        ligand = 'Inhibitor'
    if ligand not in ligand_sites:
        ligand_sites[ligand] = {}
    site = int(line.split('\t')[4])
    pdbs = line.split('\t')[7].replace('\n', '').split(';')
    if site not in ligand_sites[ligand]:
        ligand_sites[ligand][site] = pdbs
    else:
        ligand_sites[ligand][site] += pdbs
        ligand_sites[ligand][site] = list(set(ligand_sites[ligand][site]))

for ligand in ligand_sites:
    if ligand not in ligand_sites:
        print (ligand, 'does not exist in the file')
        continue
    for site in ligand_sites[ligand]:
        row = []
        row.append(ligand+' binding sites')
        row.append(int(site))
        row.append(pfam[site])
        row.append('ligand')
        row.append(ligand+'Site')
        if len(ligand_sites[ligand][site]) > 2:
            row.append(np.log2(len(ligand_sites[ligand][site])))
            row.append(np.log2(len(ligand_sites[ligand][site])))
        else:
            row.append(len(ligand_sites[ligand][site]))
            row.append(len(ligand_sites[ligand][site]))
        print (row)
        data.append(row)
# sys.exit()
'''
## Ligand binding sites
ligand_sites = {}
LIGANDS = ['ATP', 'ADP', 'MG', 'MN',
            '0WM','1LT','07J','DB8',
            '6GY','4MK','6T2','VGH',
            'P06','1N1','AQ4','E53',
            'IRE','STI','NIL','YY3',
            'LQQ','P30','BAX','B49',
            '032']
for line in open('ATP_binding_sites3.tsv', 'r'):
    if line[0] == '#':
        continue
    # Consider only the one with the specified domain
    if line.split('\t')[2] != PFAM_DOM:
        continue
    ligand = line.split('\t')[3]
    if ligand not in ligand_sites:
        ligand_sites[ligand] = {}
    if ligand not in LIGANDS:
        continue
    site = int(line.split('\t')[4])
    pdbs = line.split('\t')[7].replace('\n', '').split(';')
    if site not in ligand_sites[ligand]:
        ligand_sites[ligand][site] = pdbs
    else:
        ligand_sites[ligand][site] += pdbs
        ligand_sites[ligand][site] = list(set(ligand_sites[ligand][site]))

for ligand in LIGANDS:
    if ligand not in ligand_sites:
        print (ligand, 'does not exist in the file')
        continue
    for site in ligand_sites[ligand]:
        row = []
        row.append(ligand+'bindingSites')
        row.append(int(site))
        row.append(pfam[site])
        row.append(ligand+'site')
        row.append(ligand+'site')
        if len(ligand_sites[ligand][site]) > 2:
            row.append(np.log2(len(ligand_sites[ligand][site])))
            row.append(np.log2(len(ligand_sites[ligand][site])))
        else:
            row.append(len(ligand_sites[ligand][site]))
            row.append(len(ligand_sites[ligand][site]))
        data.append(row)
'''

## Interface binding sites
interface_sites = {}
for line in open('interface_sites.tsv', 'r'):
    if line[0] == '#':
        continue
    # Consider only the one with the specified domain
    if line.split('\t')[2] != PFAM_DOM:
        continue
    interface = line.split('\t')[3]
    site = int(line.split('\t')[4])
    pdbs = line.split('\t')[7].replace('\n', '').split(';')
    if site not in interface_sites:
        interface_sites[site] = pdbs
    else:
        interface_sites[site] += pdbs
        interface_sites[site] = list(set(interface_sites[site]))


for site in interface_sites:
    row = []
    row.append('Interface')
    row.append(int(site))
    row.append(pfam[site])
    row.append('Interface')
    row.append('Interface')
    if len(interface_sites[site]) > 2:
        row.append(np.log2(len(interface_sites[site])))
        row.append(np.log2(len(interface_sites[site])))
    else:
        row.append(len(interface_sites[site]))
        row.append(len(interface_sites[site]))
    data.append(row)

## SS of HMM
# hmm_ss = {}
for line in open('../pfam/'+PFAM_DOM+'.hmm', 'r'):
    if len(line.split()) < 20:
        continue
    if line.split()[-2] != '-' and line.split()[-3] != '-':
        continue
    ss_type = line.split()[-1].replace('\n', '')
    ss_pos = int(line.split()[0].replace('\n', ''))
    # hmm_ss[ss_pos] = ss_type
    row = []
    row.append('SS_HMM')
    row.append(ss_pos)
    row.append(pfam[ss_pos])
    row.append(ss_type)
    row.append(ss_type)
    row.append(2)
    row.append(2)
    data.append(row)

## Add PTM sites
phospho_sites = {}
acetyl_sites = {}
methyl_sites ={}
ubiq_sites = {}
sumo_sites = {}
ga_sites = {}
gl_sites = {}
for line in open('Kinase_psites4.tsv', 'r'):
    if line[0] == '#':
        continue
    # Consider only the one with the specified domain
    if line.split('\t')[2] != PFAM_DOM:
        continue
    pfam_pos = int(line.split('\t')[4])
    # num_ptmsites = int(line.split('\t')[3].replace('\n', ''))
    type_ptm = line.split('\t')[3].split('-')[1]
    if type_ptm == 'p':
        dic = phospho_sites 
    elif type_ptm == 'ac':
        dic = acetyl_sites
    elif type_ptm in ['m1', 'm2', 'm3', 'me']:
        dic = methyl_sites
    elif type_ptm in ['ub']:
        dic = ubiq_sites
    elif type_ptm in ['sm']:
        dic = sumo_sites
    elif type_ptm in ['ga']:
        dic = ga_sites
    elif type_ptm in ['gl']:
        dic = gl_sites
    else:
        print ('error:', type_ptm, 'not knonw at', pfam_pos)
        sys.exit()
    if pfam_pos not in dic:
        dic[pfam_pos] = 1
    else:
        dic[pfam_pos] += 1

count = 0
for dic, ptm_type in zip([phospho_sites, acetyl_sites, methyl_sites, ubiq_sites, sumo_sites, ga_sites, gl_sites], ['p', 'ac', 'me', 'ub', 'sm', 'ga', 'gl']):
    count += 1
    for pfam_pos in dic:
        row = []
        row.append(ptm_type+'-sites')
        row.append(int(pfam_pos))
        row.append(pfam[pfam_pos])
        row.append('PTM')
        row.append(ptm_type+'-site')
        if count in [1, 4]:
            row.append(np.log2(dic[pfam_pos])+1)
            row.append(np.log2(dic[pfam_pos])+1)
        else:
            row.append(dic[pfam_pos])
            row.append(dic[pfam_pos])
        data.append(row)

df = pd.DataFrame(data=data, columns=['Kinase', 'Pfam_Position', 'Pfam_Residue', 'Category', 'Mutation', 'Num_Samples', 'Num_Samples_Size'])
# df = df.sort_values(by=['Kinase'])
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
                "resistance": "blue",
                "ligand": "purple",
                "PTM": "grey"})
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