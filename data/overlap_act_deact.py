#!/usr/bin/env python3.10
# coding: utf-8

'''
Script to find positins where you see both activating and deactivating mutations
'''

import os, sys
sys.path.append('../ML/')
import fetchData

mydb = fetchData.connection(db_name='kinase_project2')
mydb.autocommit = True
mycursor = mydb.cursor()

mycursor.execute("SELECT pfampos, ptmtype FROM ptms\
                where pfampos!=%s", ('-',))
ptms = mycursor.fetchall()
dic_ptms = {}
for ptm_row in ptms:
    pfampos, ptmtype = ptm_row
    if ptmtype != 'p': continue
    if pfampos not in dic_ptms: dic_ptms[pfampos] = 0
    dic_ptms[pfampos] += 1

# extracts mutations
mycursor.execute("SELECT acc, wtpos, mutation, pfampos, mut_type FROM mutations\
                    where mut_type!=%s", ('neutral',))
mutations = mycursor.fetchall()
dic_hmm = {}
for mutation in mutations:
    acc, wtpos, mutation, pfampos, mut_type = mutation
    if mut_type =='activating': mut_type = 'constitutive-activation'
    if pfampos == '-': continue
    mutAA = mutation[-1]
    # print (mutAA)
    
    if pfampos not in dic_hmm: dic_hmm[pfampos] = {'constitutive-activation': 0, 'decrease': 0, 'increase': 0, 'loss': 0, 'resistance': 0}
    dic_hmm[pfampos][mut_type] += 1

for pfampos in dic_hmm:
    num_act = dic_hmm[pfampos]['constitutive-activation'] + dic_hmm[pfampos]['increase']
    num_deact = dic_hmm[pfampos]['loss'] + dic_hmm[pfampos]['decrease']
    if num_act >= 2 and num_deact >= 2 and num_act + num_deact >= 5:
        if pfampos in dic_ptms:
            print (pfampos, num_act, num_deact, dic_ptms[pfampos])
        else:
            print (pfampos, num_act, num_deact, 0)