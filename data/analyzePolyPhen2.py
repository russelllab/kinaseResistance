#!/usr/bin/env python

'''
code to prepare data to run PolyPhen2
'''

import os, sys, gzip
sys.path.insert(1, '../ML/')
import fetchData
from sklearn.metrics import matthews_corrcoef, recall_score

class Mutation:
    def __init__(self, acc, pos, wtAA, mutAA, mut_type):
        self.acc = acc
        self.pos = pos
        self.wtAA = wtAA
        self.mutAA = mutAA
        self.mut_type = mut_type
        self.prediction = None

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
    # print (acc, pos, wtAA, mutAA)
    # print (acc, mutation[0])

# get all the predictions from the PolyPhen2 output file
with gzip.open('polyphen_output.tsv.gz', 'rt') as f:
    for line in f:
        if line.startswith('#'):
            continue
        else:
            line = line.strip().split('\t')
            acc = line[0].strip()
            pos = line[1].strip()
            wtAA = line[2].strip()
            mutAA = line[3].strip()
            name = acc + '/' + wtAA+ pos + mutAA
            prediction = line[9].strip()
            # print (name, acc, pos, wtAA, mutAA, prediction)
            if name in dic_mutations:
                dic_mutations[name].prediction = prediction
            else:
                print (f'{name} not found in the DB')

y_pred = []; y_true = []
for mutation in dic_mutations:
    if dic_mutations[mutation].prediction is None: continue
    if dic_mutations[mutation].mut_type in ['resistance']: continue
    # print (mutation, dic_mutations[mutation].mut_type, dic_mutations[mutation].prediction)
    if 'damaging' in dic_mutations[mutation].prediction:
        y_pred.append(1)
    else:
        y_pred.append(0)
    if dic_mutations[mutation].mut_type in ['increase', 'activating', 'decrease', 'loss', 'resistance']:
        y_true.append(1)
    else:
        y_true.append(0)

# print (y_pred)
print (matthews_corrcoef(y_true, y_pred))
print (recall_score(y_true, y_pred))
print (recall_score(y_true, y_pred, pos_label=0))