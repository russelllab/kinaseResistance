#!/usr/bin/env python3.10
# coding: utf-8

import os, sys
sys.path.append('../ML/')
import fetchData
from tqdm import tqdm
import pickle
import pandas as pd
from cls import Kinase, Mutation
import argparse
import shutil
import gzip, sys

mydb = fetchData.connection(db_name='kinase_project2')
mydb.autocommit = True
mycursor = mydb.cursor()

mycursor.execute("SELECT acc, ptmtype, uniprotpos FROM ptms")
ptms = mycursor.fetchall()
dic_ptm = {}
for ptm_row in ptms:
    acc, ptm_type, uniprotpos = ptm_row
    uniprotpos = int(uniprotpos)
    if acc not in dic_ptm: dic_ptm[acc] = {}
    if ptm_type not in dic_ptm[acc]: dic_ptm[acc][ptm_type] = []
    dic_ptm[acc][ptm_type].append(uniprotpos)

# print (dic_ptm['P08069']['p'])
# sys.exit()

mycursor.execute("SELECT mutation, pfampos, mut_type, gene, acc, info FROM mutations\
                  WHERE pfampos!=%s", ('-',))
mutations = mycursor.fetchall()
dic_mutation_count = {}
dic_ptm_count = {}
for mutation_row in mutations:
    mutation, pfampos, mut_type, gene, acc, info = mutation_row
    if mut_type == 'neutral': continue
    # if mut_type == 'resistance': continue
    if int(pfampos) >= 795:
        uniprotpos = int(mutation[1:-1])
        if mut_type not in dic_mutation_count: dic_mutation_count[mut_type] = 0
        dic_mutation_count[mut_type] += 1
        # print (gene, mutation, pfampos, mut_type, info)
        ptm_type = '-'
        if acc in dic_ptm:
            for ptm in dic_ptm[acc]:
                if uniprotpos in dic_ptm[acc][ptm]:
                    ptm_type = ptm
                    break
        # if acc == 'P08069':
        #     print (gene, mutation, pfampos, mut_type, ptm_type)
        if ptm_type != '-':
            if ptm_type not in dic_ptm_count: dic_ptm_count[ptm_type] = 0
            dic_ptm_count[ptm_type] += 1
        # if 'C' in mutation:
        print (f'{gene}\t{mutation}\t{ptm_type}\t{pfampos}\t{mut_type}\t{info}')

print (dic_mutation_count)
print (dic_ptm_count)