#!/usr/bin/env python3.10
# coding: utf-8

import gzip
import threading
import fetchData
from tqdm import tqdm
import prepareTestData
from cls import Kinase, Mutation

AA = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

accs_to_consider = []
for line in gzip.open('all_kinases_acc.txt.gz', 'rt'):
    accs_to_consider.append(line.rstrip().lstrip())

mydb = fetchData.connection(db_name='kinase_project2')
mydb.autocommit = True
mycursor = mydb.cursor()

mycursor.execute("SELECT pfampos, name FROM positions")
hits = mycursor.fetchall()
kinases = {}
for row in hits :
    pfampos, name = row
    if pfampos == '-': continue
    # print (pfampos)
    acc = name.split('/')[0]
    if acc not in accs_to_consider: continue
    if acc not in kinases: kinases[acc] = []
    kinases[acc].append(name.split('/')[1])

def call_activark(acc):
    outputFile = 'outputs/'+acc+'.txt'
    prepareTestData.predict(5000, 'inputs/'+acc+'.txt.gz', outputFile = outputFile)

count = 0
for acc in tqdm(kinases):
    count += 1
    mutations = ''
    for num, case in enumerate(kinases[acc], start=1):
        wtAA = case[0]
        position = int(case[1:])
        for mutAA in AA:
            if mutAA == wtAA: continue
            # if mutAA != 'A': continue
            mutations += acc + '/' + wtAA + str(position) + mutAA + '\n'
        #if num % 1 == 0:
        #    break
    
    gzip.open('inputs/'+acc+'.txt.gz', 'wt').write(mutations)
    call_activark(acc)
    # thread = threading.Thread(target=call_activark, args=(acc,))
    # thread.start()
    # prepareTestData.predict('inputs/'+acc+'.txt')
    # if count == 5: break
    # while threading.active_count() > 25:
    #     pass
    
mydb.close()
