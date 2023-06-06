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
sys.path.insert(1, '../ML/')
import fetchData

import numpy as np
from matplotlib_venn import venn3, venn3_circles

mydb = fetchData.connection(db_name='kinase_project2')
mydb.autocommit = True
mycursor = mydb.cursor()

mycursor.execute("SELECT mutation, pfampos, mut_type FROM mutations\
		 		where pfampos!=%s and pfampos!=%s", ('-', 'neutral'))
mutations = mycursor.fetchall()
a =[]
d=[]
r=[]
for mut_row in mutations:
	mutation, pfampos, mut_type = mut_row
	if mut_type in ['loss', 'decrease']:
		d.append(pfampos)
		# d.append(mutation)
	elif mut_type in ['activating', 'increase']:
		a.append(pfampos)
		# a.append(mutation)
	elif mut_type in ['resistance']:
		r.append(pfampos)
		# r.append(mutation)

print (len(a), len(d), len(r))
print (len(list(set(a+d+r))))

# Make a Basic Venn
# a = [1,1,1]
# b = [1,2,3,2,1]
# c = [1,2,1,2,1,2,3,4]
# v = venn3(subsets=(1, 2, 2, 2, 2, 2, 2), set_labels = ('A', 'B', 'C'))
v = venn3([set(a), set(d), set(r)], set_labels = ('Activating', 'Deactivating', 'Resistance'))
 
# Custom it
# v.get_patch_by_id('100').set_alpha(1.0)
# v.get_patch_by_id('100').set_color('white')
# v.get_label_by_id('100').set_text('Unknown')
# v.get_label_by_id('A').set_text('Set "A"')
# c = venn3_circles(subsets=(1, 1, 1, 1, 1, 1, 1), linestyle='dashed')
c = venn3_circles([set(a), set(d), set(r)], linestyle='dashed')
# c[0].set_lw(1.0)
# c[0].set_ls('dotted')
 
# Add title and annotation
plt.title("Distribution of mutated positions")
# plt.annotate('Unknown set', xy=v.get_label_by_id('100').get_position() - np.array([0, 0.05]), xytext=(-70,-70),
# ha='center', textcoords='offset points', bbox=dict(boxstyle='round,pad=0.5', fc='gray', alpha=0.1),
# arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.5',color='gray'))
 
# Show it
# plt.show()
# Save it
plt.savefig('pfam_venn.png', format='png', dpi=300, bbox_inches='tight')
plt.savefig('pfam_venn.svg', format='png', dpi=300, bbox_inches='tight')