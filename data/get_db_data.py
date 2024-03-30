#!/usr/bin/env python3
# coding: utf-8

'''
A script to cluster kinase resistance, activating,
and deactivating mutations
'''
import os, sys, gzip
from turtle import position
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import plotly.graph_objects as go
from tqdm import tqdm
sys.path.insert(1, '../ML/')
import fetchData

import numpy as np
from matplotlib_venn import venn3, venn3_circles

mydb = fetchData.connection(db_name='kinase_project2')
mydb.autocommit = True
mycursor = mydb.cursor()

mycursor.execute("SELECT table_name FROM information_schema.tables WHERE table_schema = 'public'")
tables = mycursor.fetchall()
print (tables)
for table in tqdm(tables):
    text = ''
    table = table[0]
    # print (table)
    # mycursor.execute("Select * FROM people LIMIT 0")
    mycursor.execute("SELECT * FROM "+table+" LIMIT 0")
    data = mycursor.description
    values = []
    for h in data:
        # print (h[0])
        values.append(str(h[0]))
    text += '\t'.join(values) + '\n'
    mycursor.execute("SELECT * FROM "+table)
    data = mycursor.fetchall()
    # print (data)
    for row in data:
        values = []
        for value in row:
            values.append(str(value))
        text += '\t'.join(values) + '\n'
        # print (text)
    with open('dataTables/'+table+'.txt.gz', 'wt') as f:
        f.write(text)
    # sys.exit()

print ('done')