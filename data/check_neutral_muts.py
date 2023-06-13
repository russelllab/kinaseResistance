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

class Kinase:
    def __init__(self, acc, gene):
        self.acc = acc
        self.gene = gene
        self.fasta = ''
        self.mutations = []

kinases = {}
## read neutral mutations
for line in open('20230525_AllNeutrals_1pSNP_5pHom.txt', 'r'):
    if line.startswith('Uniprot'): continue
    acc, gene, mutation = line.strip().split()[0], line.strip().split()[1], line.strip().split()[2]
    if acc not in kinases:
        kinases[acc] = Kinase(acc, gene)
    kinases[acc].mutations.append(mutation)

# print (kinases.keys())

mydb = fetchData.connection(db_name='kinase_project2')
mydb.autocommit = True
mycursor = mydb.cursor()

mycursor.execute("SELECT acc, gene, fasta from kinases")

hits = mycursor.fetchall()
for hit in hits:
    acc, gene, fasta = hit
    # print (acc, gene, fasta)
    if acc not in kinases: continue
    kinases[acc].fasta = fasta
    # break

for acc in kinases:
    gene = kinases[acc].gene
    for mutation in kinases[acc].mutations:
        pos = mutation[1:-1]
        wt = mutation[0]
        mut = mutation[-1]
        if pos.isnumeric() == False: continue
        if int(pos) > len(kinases[acc].fasta):
            print (f'{acc}\t{gene}\t{mutation}\tProtein length: {len(kinases[acc].fasta)}')
        elif kinases[acc].fasta[int(pos)-1] != wt:
            print (f'{acc}\t{gene}\t{mutation}\t{kinases[acc].fasta[int(pos)-1]}')