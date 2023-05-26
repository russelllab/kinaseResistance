#!/usr/bin/env python3
# coding: utf-8

import sys
sys.path.append('../ML/')
import fetchData

'''
To make comparisons between the different datasets.
'''

dic1 = {}
dic2 = {}
acc2gene = {}
for dic, db_name in zip([dic1, dic2], ['kinase_project', 'kinase_project2']):
    mydb = fetchData.connection(db_name=db_name)
    mydb.autocommit = True
    mycursor = mydb.cursor()

    # select all mutations from the database
    mycursor.execute("SELECT mutation, mut_type, pfampos, acc, gene FROM mutations")
    hits = mycursor.fetchall()
    for hit in hits:
        mutation, mut_type, pfampos, acc, gene = hit
        if acc not in acc2gene: acc2gene[acc] = gene
        if pfampos == '-': continue
        if acc not in dic:
            dic[acc] = {}
        if mutation not in dic[acc]:
            dic[acc][mutation] = []
        if mut_type not in dic[acc][mutation]:
            dic[acc][mutation].append(mut_type)

for acc in dic2:
    if acc not in dic1:
        print (f'{acc} {acc2gene[acc]} not in kinase_project. {dic2[acc]}')
        continue
    for mutation in dic2[acc]:
        if mutation not in dic1[acc]:
            if len(dic2[acc][mutation])>1 or 'N' not in dic2[acc][mutation]:
                print (f'{acc}/{mutation} {acc2gene[acc]} {dic2[acc][mutation]} not in {acc} of kinase_project')
            continue
        if dic2[acc][mutation] != dic1[acc][mutation]:
            print (f'{acc} {acc2gene[acc]} {mutation} old_data:{dic1[acc][mutation]} new_data:{dic2[acc][mutation]}')
            continue
