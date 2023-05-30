#!/usr/bin/env python

'''
code to prepare data to run PolyPhen2
'''

import os, sys
sys.path.insert(1, '../ML/')
import fetchData

mydb = fetchData.connection(db_name='kinase_project2')
mydb.autocommit = True
mycursor = mydb.cursor()

# get all the variants
mycursor.execute("SELECT mutation, mut_type, acc, gene FROM mutations where\
                    pfampos!=%s", ('-',))
mutations = mycursor.fetchall()
for mutation in mutations:
    # print (mutation)
    wtAA = mutation[0][0]
    mutAA = mutation[0][-1]
    pos = mutation[0][1:-1]
    acc = mutation[2]
    print (acc, pos, wtAA, mutAA)
    # print (acc, mutation[0])