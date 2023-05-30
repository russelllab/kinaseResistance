#!/usr/bin/env python

'''
code to analyze to SIFT4G
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

# get all the predictions from the SIFT4G output file
mechXpath = '/net/home.isilon/ds-russell/mechismoX/analysis/mutations/data/VLatest/'
sift_dic = {}
for name in dic_mutations:
    acc = name.split('/')[0]
    if acc in sift_dic: continue
    sift_dic[acc] = {}
    for line in gzip.open(mechXpath+acc[:4]+acc+'P15056.annotations.txt.gz', 'rt'):
        if line.startswith('UniProtID'): continue
        mutation = line.strip().split('\t')[0]+'/'+line.split('\t')[1]+line.split('\t')[2]+line.split('\t')[3]
        siftPred = 'damaging' if 'True' in line.strip().split('\t')[19] else 'benign'
        sift_dic[acc][mutation] = siftPred
        