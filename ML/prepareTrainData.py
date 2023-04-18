#!/usr/bin/env python
# coding: utf-8

## Develop training set data
import numpy as np
import scipy as sp
import os, sys, gzip
# import seaborn as sns
# import pandas as pd
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
import matplotlib.pyplot as plt
from sklearn import decomposition
from sklearn.preprocessing import MinMaxScaler
tf.random.set_seed(69)
from cls import Kinase, Mutation
import fetchData

PTM_TYPES = ['ac', 'gl', 'm1', 'm2', 'm3', 'me', 'p', 'sm', 'ub']
AA = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
WS = 5

exceptions= ['Q9Y2K2', 'Q15303', 'Q9UIK4', 'P33981', 'P35916',
             'Q99683', 'Q8IVW4', 'Q9H792', 'Q9P286', 'Q86Z02',
             'Q8TF76', 'Q96L34', 'Q13308', 'Q9UK32', 'Q15772',
             'P51617', 'Q9Y3S1', 'Q9C098', 'Q6VAB6', 'P21127',
             'Q13557', 'Q6ZMQ8', 'Q6P0Q8', 'Q8IZE3', 'P51957',
             'O60229', 'Q96RG2', 'Q5VST9', 'Q8WZ42', 'O75962',
             'O95835', 'Q13535']

kinases = {}

mydb = fetchData.connection()
mydb.autocommit = True
mycursor = mydb.cursor()
mycursor.execute("select fasta from kinases where acc=%s", ('P06493',))
fasta = mycursor.fetchone()[0]
print (fasta)

# fetchData.fetchFasta(kinases, Kinase)
fetchData.fetchFasta(kinases, Kinase, mycursor)
fetchData.fetchGroup(kinases, Kinase)
# print (kinases['P00533'].group)
hmmPkinase = fetchData.fetchPkinaseHMM(mycursor) # hmmPosition > AA > bit-score
# print (hmmPkinase[30]['K'])
fetchData.fetchHmmsearch(kinases, Kinase)
# fetchData.dsspScores(kinases, Kinase)
# # print (kinases['Q9NYV4'].burr[3])
# # print (kinases['Q92772'].dihedral)
# fetchData.iupredScores(kinases, Kinase)
# fetchData.homologyScores(kinases, Kinase)
# sys.exit()

# #print (kinases['Q9NYV4'].mechismo)
# data = []
# for acc in kinases:
#     #print (kinases[acc].domains)
#     for domainNum in kinases[acc].domains:
#         #print (domainNum)
#         data.append(len(kinases[acc].domains[domainNum]))

# #print (data)
# df = pd.DataFrame(data = data, columns=['Length'])
# sns.histplot(data=df, x="Length")

'''
def fetchPkinase(acc, domainNum):
    # A function to take acc and alignment cutoff values to return
    # the HMM bitscores
    
    row = []
    for hmmPosition in range(1,265):
        if hmmPosition in kinases[acc].domains[domainNum]:
            SeqPosition = kinases[acc].domains[domainNum][hmmPosition]
            aa = kinases[acc].fasta[SeqPosition-1]
            value = float(hmmPkinase[hmmPosition][aa])
        else:
            value = 3

        row.append(value)
    #print (len(row))
    kinases[acc].hmm[domainNum] = row
    return row
'''

def fetchStrucFeat(acc, domainNum):
    data = []
    df = pd.DataFrame()
    for dic, name in zip([
                        kinases[acc].dihedral,
                        kinases[acc].sec,
                        kinases[acc].burr,
                        kinases[acc].access,
                        kinases[acc].iupred,
                        kinases[acc].mechismo
                        ],[
                        'dihedral',
                        'access',
                        'burr',
                        'sec',
                        'iupred',
                        'mechismo'
                        ]):
        row = []
        print (kinases[acc].dihedral)
        sys.exit()
        for hmmPosition in range(1,265):
            if hmmPosition in kinases[acc].domains[domainNum]:
                SeqPosition = kinases[acc].domains[domainNum][hmmPosition]
                residue = kinases[acc].fasta[SeqPosition-1]
                #print (SeqPosition)
                try:
                    value = dic[SeqPosition][residue]
                except:
                    print (acc, SeqPosition, len(dic), residue)
                    sys.exit()
            else:
                value = 3
            row.append(value)
        #f = pd.DataFrame(data, columns=AA)
        df[name] = row
    #print (df)

    return df

'''Map sequence to Pfam for all canonical kinases'''
seq2pfam = {}
mycursor.execute("select acc, uniprotpos, pfampos from positions")
for row in mycursor.fetchall():
    acc = row[0]
    seqPos = str(row[1])
    pfamPos = str(row[2])
    if acc not in seq2pfam: seq2pfam[acc] = {}
    seq2pfam[acc][seqPos] = pfamPos
'''
for line in gzip.open('../data/humanKinasesHitsHmmsearchMappings.tsv.gz', 'rt'):
    if line[0] =='#':
        continue
    acc = line.split('\t')[0].split('|')[1]
    seqPos = str(line.split('\t')[2])
    pfamPos = str(line.replace('\n', '').split('\t')[4])
    if acc not in seq2pfam: seq2pfam[acc] = {}
    seq2pfam[acc][seqPos] = pfamPos
'''

'''Fetch test mutation data'''
for line in open('test_mutations.txt', 'r'):
    gene = line.split()[1]
    acc = line.split()[0]
    mutation = line.split()[2]
    wtAA = mutation[0]
    mutAA = mutation[-1]
    position = str(mutation[1:-1])
    mut_type = line.split()[3].replace('\n', '')
    dataset = 'test'
    # print (acc, kinases[acc].gene, wtAA, position, mutAA)
    if position in seq2pfam[acc]:
        kinases[acc].mutations[mutation] = Mutation(mutation, mut_type, acc, dataset)
        kinases[acc].mutations[mutation].positionHmm = seq2pfam[acc][position]
    # pkinase_act_deact_res[mut_type].append(kinases[acc].mutations[mutation].positionHmm)

pkinase_act_deact_res = {'A': [], 'D': [], 'R': [], 'N': []}
'''Fetch act/deact mutation data'''
mycursor.execute("select acc, gene, mutation, wtaa, mutaa, wtpos, mut_type from mutations")
for row in mycursor.fetchall():
    acc = row[0]
    gene = row[1]
    mutation = row[2]
    wtAA = row[3]
    mutAA = row[4]
    position = str(row[5])
    mut_type = row[6]
    dataset = 'train'
    # print (acc, kinases[acc].gene, wtAA, position, mutAA)
    if position in seq2pfam[acc]:
        kinases[acc].mutations[mutation] = Mutation(mutation, mut_type, acc, dataset)
        kinases[acc].mutations[mutation].positionHmm = seq2pfam[acc][position]
        pkinase_act_deact_res[mut_type].append(kinases[acc].mutations[mutation].positionHmm)
'''
for line in open('../AK_mut_w_sc_feb2023/act_deact_v2.tsv', 'r'):
    if line.split()[0] == 'uniprot_name': continue
    gene = line.split('\t')[0]
    acc = line.split('\t')[1]
    wtAA = line.split('\t')[2].replace(',', '')
    mutAA = line.split('\t')[4].replace(',', '')
    if len(wtAA) > 1 or len(mutAA) > 1: continue
    position = str(line.split('\t')[3])
    mut_type = line.split('\t')[5]
    # print (acc, kinases[acc].gene, wtAA, position, mutAA)
    mutation = wtAA + position + mutAA
    if acc not in seq2pfam:
        continue
    if mutation not in kinases[acc].mutations:
        dataset = 'train'
        kinases[acc].mutations[mutation] = Mutation(mutation, mut_type, acc, dataset)
        # print (acc, gene, position, wtAA, mutAA)
        kinases[acc].mutations[mutation].positionHmm = seq2pfam[acc][position]
        pkinase_act_deact_res[mut_type].append(kinases[acc].mutations[mutation].positionHmm)
'''

'''Fetch resistant mutation data'''
'''
for line in gzip.open('../KA/resistant_mutations_Mar_2023.tsv.gz', 'rt'):
    if line[0] == '#': continue
    actual_gene = line.split('\t')[0]
    if '_' in actual_gene: gene = actual_gene.split('_')[0]
    else: gene = actual_gene
    acc = line.split('\t')[2]
    cosmic_mutation = line.split('\t')[1]
    wtAA = cosmic_mutation[0]
    mutAA = cosmic_mutation[-1]
    if mutAA == 'X': continue
    if len(wtAA) > 1 or len(mutAA) > 1: continue
    uniprot_position = line.split('\t')[5].replace('\n', '')
    mutation = wtAA + uniprot_position + mutAA
    mut_type = 'R'
    # print (acc, kinases[acc].gene, wtAA, position, mutAA)
    if mutation not in kinases[acc].mutations:
        dataset = 'train'
        kinases[acc].mutations[mutation] = Mutation(mutation, mut_type, acc, dataset)
        try:
            kinases[acc].mutations[mutation].positionHmm = seq2pfam[acc][uniprot_position]
            pkinase_act_deact_res[mut_type].append(kinases[acc].mutations[mutation].positionHmm)
        except:
            print (acc, actual_gene, uniprot_position, mutation, cosmic_mutation)
            sys.exit()
    else: kinases[acc].mutations[mutation].mut_types.append(mut_type)
    # kinases[acc].mutations[wtAA+position+mutAA] = Mutation(wtAA+position+mutAA, mut_type)
'''

'''Fetch resistant mutation data'''
for line in open('../AK_mut_w_sc_feb2023/res_mut_v3_only_subs_KD_neighb.tsv', 'r'):
    if line.split('\t')[0] == 'uniprot_id': continue
    # actual_gene = line.split('\t')[0]
    # if '_' in actual_gene: gene = actual_gene.split('_')[0]
    # else: gene = actual_gene
    acc = line.split('\t')[0]
    # cosmic_mutation = line.split('\t')[1]
    wtAA = line.split('\t')[1]
    mutAA = line.split('\t')[3]
    if mutAA == 'X': continue
    if len(wtAA) > 1 or len(mutAA) > 1: continue
    uniprot_position = line.split('\t')[2].replace('\n', '')
    mutation = wtAA + uniprot_position + mutAA
    mut_type = 'R'
    if acc not in seq2pfam:
        continue
    if uniprot_position not in seq2pfam[acc]:
        print (f'{uniprot_position} seems to be outside the domain in {acc} and reported {mut_type}')
        continue
    # print (acc, kinases[acc].gene, wtAA, position, mutAA)
    if mutation not in kinases[acc].mutations:
        dataset = 'train'
        kinases[acc].mutations[mutation] = Mutation(mutation, mut_type, acc, dataset)
        try:
            kinases[acc].mutations[mutation].positionHmm = seq2pfam[acc][uniprot_position]
            pkinase_act_deact_res[mut_type].append(kinases[acc].mutations[mutation].positionHmm)
        except:
            print (acc, uniprot_position, mutation, mut_type)
            sys.exit()
    else: kinases[acc].mutations[mutation].mut_types.append(mut_type)
    # kinases[acc].mutations[wtAA+position+mutAA] = Mutation(wtAA+position+mutAA, mut_type)

'''Fetch neutral mutation data'''
for line in open('../AK_mut_w_sc_feb2023/nat_mut_tidy_v2_march2023.tsv', 'r'):
    if line.split('\t')[1] == 'UniProtID': continue
    acc = line.split('\t')[1]
    wtAA = line.split('\t')[2]
    mutAA = line.split('\t')[4]
    if mutAA == 'X': continue
    if len(wtAA) > 1 or len(mutAA) > 1: continue
    uniprot_position = line.split('\t')[3].replace('\n', '')
    mutation = wtAA + uniprot_position + mutAA
    mut_type = 'N'
    if acc not in seq2pfam:
        continue
    if uniprot_position not in seq2pfam[acc]:
        print (f'{uniprot_position} seems to be outside the domain and reported {mut_type}')
        print (seq2pfam[acc])
        continue
    # print (acc, kinases[acc].gene, wtAA, position, mutAA)
    if mutation not in kinases[acc].mutations:
        dataset = 'train'
        kinases[acc].mutations[mutation] = Mutation(mutation, mut_type, acc, dataset)
        try:
            kinases[acc].mutations[mutation].positionHmm = seq2pfam[acc][uniprot_position]
            pkinase_act_deact_res[mut_type].append(kinases[acc].mutations[mutation].positionHmm)
        except:
            print (acc, uniprot_position, mutation, mut_type)
            sys.exit()
    else: kinases[acc].mutations[mutation].mut_types.append(mut_type)
    # kinases[acc].mutations[wtAA+position+mutAA] = Mutation(wtAA+position+mutAA, mut_type)

'''Fetch PTM data'''
hmmPTM = {}
mycursor.execute("select acc, ptmtype, uniprotpos, pfampos from ptms")
for row in mycursor.fetchall():
    acc = row[0]
    ptm_type = row[1]
    uniprot_position = row[2]
    hmm_position = int(row[3])
    if ptm_type not in kinases[acc].ptm:
        kinases[acc].ptm[ptm_type] = []
    kinases[acc].ptm[ptm_type].append(int(uniprot_position))
    if hmm_position not in hmmPTM:
        hmmPTM[hmm_position] = []
    hmmPTM[hmm_position].append(ptm_type)
# print (hmmPTM[200])
       
'''
# for line in open('../data/Kinase_psites4.tsv', 'r'):
# for line in open('../data/Kinase_psites_trimmed.tsv', 'r'):
for line in open('../data/Kinase_psites_hits_split_trimmed.tsv', 'r'):
    if line[0] == '#': continue
    if line.split()[2] != 'humanKinasesHitsSplitTrimmed': continue
    acc = line.split('\t')[0]
    if acc not in kinases: continue
    position = int((line.split('\t')[3].split('-')[0])[1:])
    ptm_type = line.split('\t')[3].split('-')[1]
    hmm_position = int(line.split('\t')[4])
    if ptm_type not in kinases[acc].ptm:
        kinases[acc].ptm[ptm_type] = []
    kinases[acc].ptm[ptm_type].append(position)
    if hmm_position not in hmmPTM:
        hmmPTM[hmm_position] = []
    hmmPTM[hmm_position].append(ptm_type)
    # print (acc, kinases[acc].gene, position, ptm_type)
# print (hmmPTM[141])
'''

'''Make training matrix'''
trainMat = 'Acc\tGene\tMutation\tDataset\t'
trainMat += 'hmmPos\thmmSS\thmmScoreWT\thmmScoreMUT\thmmScoreDiff\t'
trainMat += 'Phosphomimic\t'
trainMat += 'ChargesWT\tChargesMUT\tChargesDiff\t'
startWS = int((WS-1)/2) * -1
endWS = int((WS-1)/2)
for position in range(startWS, endWS+1):
    trainMat += ('_'+str(position)+'\t').join(PTM_TYPES) + '_'+str(position)+'\t'
    trainMat += ('_'+str(position)+'_pfam\t').join(PTM_TYPES) + '_' + str(position) + '_pfam\t'
trainMat += '_WT\t'.join(AA) + '_WT\t'
trainMat += '_MUT\t'.join(AA) + '_MUT\t'
# trainMat += '\t'.join(['allHomologs','exclParalogs','specParalogs','orthologs','bpso','bpsh']) + '\t'
for position in range(startWS, endWS+1):
    trainMat += ('_'+str(position)+'\t').join(['A', 'D', 'R']) + '_'+str(position)+'\t'
# trainMat += '_known\t'.join(['A', 'D', 'R']) + '_known\t'
trainMat += 'MUT_TYPE\n'
# print (trainMat)
# print ('_WT\t'.join(AA) + '\t')
# print (AA)
# sys.exit()
data = []
mut_types_colors = []
for acc in kinases:
    if len(kinases[acc].mutations) == 0:
        continue
    for mutation in kinases[acc].mutations:
        row = []
        mutation_obj = kinases[acc].mutations[mutation]
        position = mutation_obj.position
        mutAA = mutation_obj.mutAA
        wtAA = mutation_obj.wtAA
        mut_types = ''.join(np.sort(list(set(mutation_obj.mut_types))))
        hmmPos, hmmScoreWT, hmmScoreMUT, hmmSS = fetchData.getHmmPkinaseScore(acc, wtAA, position, mutAA, kinases, hmmPkinase)
        ptm_row = fetchData.getPTMscore(acc, position, kinases, hmmPTM, WS)
        aa_row = fetchData.getAAvector(wtAA, mutAA)
        # homology_row = fetchData.getHomologyScores(acc, wtAA, position, mutAA, kinases)
        is_phosphomimic = kinases[acc].mutations[mutation].checkPhosphomimic()
        charges_row = kinases[acc].mutations[mutation].findChangeInCharge()
        adr_row = fetchData.getADRvector(acc, position, kinases, pkinase_act_deact_res, WS)
        print (
            acc +'\t'+ mutation +'\t'+ str(hmmPos) +'\t'+
            str(hmmScoreWT)+'\t' +str(hmmScoreMUT)+'\t'+ ','.join(ptm_row) + '\t' +
            ','.join(aa_row) + '\t' + '\t'.join(mut_types)
            )
        # row.append(mutation_obj.dataset)
        # Prepare rows for the numpy data
        row.append(int(hmmPos))
        row.append(str(hmmSS))
        row.append(float(hmmScoreWT))
        row.append(float(hmmScoreMUT))
        row.append(float(hmmScoreMUT)-float(hmmScoreWT))
        row.append(is_phosphomimic)
        row += [int(item) for item in charges_row]
        row += [int(item) for item in ptm_row]
        row += [int(item) for item in aa_row]
        # row += [int(item) for item in homology_row]
        row += [int(item) for item in adr_row]
        data.append(row)

        # Prepare rows for writing the numpy data
        trainMat += acc + '\t' + kinases[acc].gene + '\t' + mutation + '\t' + mutation_obj.dataset + '\t'
        trainMat += str(hmmPos) + '\t' + str(hmmSS) + '\t' + str(hmmScoreWT) + '\t' + str(hmmScoreMUT) + '\t' + str(hmmScoreMUT-hmmScoreWT) + '\t'
        trainMat += str(is_phosphomimic) + '\t'
        trainMat += '\t'.join([str(item) for item in charges_row]) + '\t'
        trainMat += '\t'.join([str(item) for item in ptm_row]) + '\t'
        trainMat += '\t'.join([str(item) for item in aa_row]) + '\t'
        # trainMat += '\t'.join([str(item) for item in homology_row]) + '\t'
        trainMat += '\t'.join([str(item) for item in adr_row]) + '\t'
        trainMat += mut_types + '\n'

        if mut_types == 'A':
            mut_types_colors.append('green')
        elif mut_types == 'D':
            mut_types_colors.append('red')
        elif mut_types == 'N':
            mut_types_colors.append('cyan')
        elif mut_types == 'R':
            mut_types_colors.append('blue')
        else:
            mut_types_colors.append('violet')

gzip.open('trainDataFromHitsSplitTrimmedAln.tsv.gz', 'wt').write(trainMat)
data = np.array(data)
scaler = MinMaxScaler()
scaler.fit(data)
data = scaler.transform(data)
# print (trainMat)
# sys.exit()

pca = decomposition.PCA(n_components=2)
pca.fit(data)
X = pca.transform(data)

fig = plt.figure(1, figsize=(4, 3))
plt.clf()

ax = fig.add_subplot(111)
# ax = fig.add_subplot(111, projection="3d", elev=48, azim=134)
ax.set_position([0.1, 0.1, 0.8, 0.8])

plt.cla()

# ax.scatter(X[:, 0], X[:, 1], X[:, 2], c=mut_type_colors, cmap=plt.cm.nipy_spectral, edgecolor="k")
ax.scatter(X[:, 0], X[:, 1], c=mut_types_colors, cmap=plt.cm.nipy_spectral, edgecolor="k")
plt.show()
# plt.savefig('pca_plot.png')
