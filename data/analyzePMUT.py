#!/usr/bin/env python

'''
code to analyze to SIFT4G
'''

import os, sys, gzip
sys.path.insert(1, '../ML/')
import fetchData
from sklearn.metrics import matthews_corrcoef, recall_score, roc_auc_score

class Mutation:
    def __init__(self, acc, pos, wtAA, mutAA, mut_type):
        self.acc = acc
        self.pos = pos
        self.wtAA = wtAA
        self.mutAA = mutAA
        self.mut_type = mut_type
        self.prediction = None
        self.prob = None
    def show(self):
        print (self.acc, self.pos, self.wtAA, self.mutAA, self.mut_type, self.prediction, self.prob)

mydb = fetchData.connection(db_name='kinase_project2')
mydb.autocommit = True
mycursor = mydb.cursor()

# get all the variants from the DB
mycursor.execute("SELECT mutation, mut_type, acc, gene FROM mutations where\
                    pfampos!=%s", ('-'))
mutations = mycursor.fetchall()
dic_mutations = {}
for mutation in mutations:
    # print (mutation)
    wtAA = mutation[0][0]
    mutAA = mutation[0][-1]
    pos = mutation[0][1:-1]
    acc = mutation[2]
    mut_type = mutation[1]
    name = acc + '/' + mutation[0]
    if name not in dic_mutations:
        dic_mutations[name] = Mutation(acc, pos, wtAA, mutAA, mut_type)
    else:
        dic_mutations[name].mut_type += mut_type
    # print (acc, pos, wtAA, mutAA)
    # print (acc, mutation[0])

# get all the predictions from the SIFT4G output file
mechXpath = '/net/home.isilon/ds-russell/mechismoX/analysis/mutations/data/VLatest/'
pmut_dic = {}
pmut_output = ''
for name in dic_mutations:
    acc = name.split('/')[0]
    if acc in pmut_dic: continue
    pmut_dic[acc] = {}
    for line in gzip.open(mechXpath+acc[:4]+'/'+acc+'.annotations.txt.gz', 'rt'):
        if line.startswith('UniProtID'): continue
        mutation = line.split('\t')[1].strip()+line.split('\t')[2].strip()+line.split('\t')[3].strip()
        if line.split('\t')[19].lstrip().rstrip() in ['-', '()']: continue
        pmutPred = 'damaging' if 'True' in line.strip().split('\t')[19] else 'benign'
        # print (acc, mutation, line.split('\t')[19])
        pmut_output += acc + '\t' + mutation + '\t' + line.split('\t')[19].lstrip().rstrip() + '\n'
        pmutProb = float(line.strip().split('\t')[19].split('(')[1].split(')')[0])
        pmut_dic[acc][mutation] = pmutPred
        #print (acc+'/'+mutation)
        if acc+'/'+mutation in dic_mutations:
                dic_mutations[acc+'/'+mutation].prediction = pmutPred
                dic_mutations[acc+'/'+mutation].prob = pmutProb

gzip.open('pmut_output.tsv.gz', 'wt').write(pmut_output)

MUT_TYPES = {
            'A': ['activating', 'increase'],
            'D': ['decrease', 'loss'],
            'R': ['resistance'],
            'T': ['activatingresistance'],
            'N': ['neutral']
            }
text = ''
for mut_type in MUT_TYPES:
    names = []
    y_pred = []; y_true = []; y_prob=[]; y_act_deact_or_neutral = []
    for mutation in dic_mutations:
        if dic_mutations[mutation].prediction is None: continue
        kinase_mut_type = dic_mutations[mutation].mut_type
        if kinase_mut_type not in MUT_TYPES[mut_type]: continue
        # print (mutation, dic_mutations[mutation].mut_type, dic_mutations[mutation].prediction)
        if 'damaging' in dic_mutations[mutation].prediction:
            y_pred.append(1)
        else:
            y_pred.append(0)
        y_prob.append(float(dic_mutations[mutation].prob))
        if kinase_mut_type in MUT_TYPES['N']:
            y_act_deact_or_neutral.append(0)
        else:
            y_act_deact_or_neutral.append(1)
        if kinase_mut_type in MUT_TYPES[mut_type]:
            y_true.append(1)
        else:
            y_true.append(0)
        names.append(mutation)
    # print (f'MCC: {matthews_corrcoef(y_true, y_pred)}')
    print (mut_type, len(y_pred), f'REC: {recall_score(y_true, y_pred)}')
    # print (f'AUC: {roc_auc_score(y_true, y_pred)}')
    if mut_type == 'T':
        print (y_true)
        print (y_pred)

    for name, y_p, y_t in zip(names, y_prob, y_act_deact_or_neutral):
        if 'V600E' in name: dic_mutations[name].show()
        text += str(name) + '\t' + str(dic_mutations[name].mut_type) + '\t' + str(y_p) + '\t' + str(y_t) + '\n'

open('pmut_roc.txt', 'w').write(text)
# print (y_pred)
# print (f'MCC: {matthews_corrcoef(y_true, y_pred)}')
# print (f'REC: {recall_score(y_true, y_pred)}')
# print (f'SPE: {recall_score(y_true, y_pred, pos_label=0)}')
# print (f'AUC: {roc_auc_score(y_true, y_pred)}')
