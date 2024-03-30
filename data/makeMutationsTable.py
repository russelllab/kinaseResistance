#!/usr/bin/env python3
# coding: utf-8

'''
A script to make tables of mutations
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

class HMM:
	def __init__(self, pfampos) -> None:
		self.pfampos = pfampos
		self.mut_types = {}
		self.ptm_types = 0
		self.bitscore = None

mydb = fetchData.connection(db_name='kinase_project2')
mydb.autocommit = True
mycursor = mydb.cursor()
'''
# make mutations table
mycursor.execute("SELECT * FROM mutations")
mutations = mycursor.fetchall()
text = 'UniProtAccession\tGeneName\tMutation\tWT\tPosition\tMUT\tAlignmentPosition\tMutationType\tDescription\tPubMedID\tSource\n'
for mut_row in mutations:
    # print (len(mut_row))
    # print (mut_row)
    mutation_id, mutation, wtaa, wtpos, mutaa, pfampos, mut_type, acc, gene, info, pubmed, source = mut_row
    text += acc + '\t'
    text += gene + '\t'
    text += mutation + '\t'
    text += wtaa + '\t'
    text += str(wtpos) + '\t'
    text += mutaa + '\t'
    text += str(pfampos) + '\t'
    text += mut_type + '\t'
    text += info + '\t'
    text += pubmed.replace(',',';') + '\t'
    text += source + '\n'

open('mutations_table.tsv', 'w').write(text)
sys.exit()
# make ptms table
mycursor.execute("SELECT * FROM ptms")
ptms = mycursor.fetchall()
text = 'UniProtAccession\tGeneName\tPosition\tAlignmentPosition\tPTM\n'
for ptm_row in ptms:
    # print (len(ptm_row))
    wtpos, wtaa, pfampos,pfamaa, acc, gene, ptm_type, name = ptm_row
    text += acc + '\t'
    text += gene + '\t'
    text += str(wtpos) + '\t'
    text += str(pfampos) + '\t'
    text += ptm_type + '\n'

open('ptms_table.tsv', 'w').write(text)

# make HMM table
mycursor.execute("SELECT * FROM hmm")
hmm = mycursor.fetchall()
text = 'hmmPosition\thmmAminoAcid\talignmentPosition\t'
text += 'EmissionScore_A\tEmissionScore_C\tEmissionScore_D\t'
text += 'EmissionScore_E\tEmissionScore_F\tEmissionScore_G\t'
text += 'EmissionScore_H\tEmissionScore_I\tEmissionScore_K\t'
text += 'EmissionScore_L\tEmissionScore_M\tEmissionScore_N\t'
text += 'EmissionScore_P\tEmissionScore_Q\tEmissionScore_R\t'
text += 'EmissionScore_S\tEmissionScore_T\tEmissionScore_V\t'
text += 'EmissionScore_W\tEmissionScore_Y\n'
for hmm_row in hmm:
    hmmpos, hmm_aa, alignpos = hmm_row[0], hmm_row[1], hmm_row[3]
    if hmmpos == '-': continue
    emm_probs = hmm_row[4:]
    text += str(hmmpos) + '\t'
    text += hmm_aa + '\t'
    text += str(alignpos) + '\t'
    for emm_prob in emm_probs:
        text += str(emm_prob) + '\t'
    text += '\n'

open('hmm_table.tsv', 'w').write(text)

# make alignment table
mycursor.execute("SELECT * FROM alignment")
alignments = mycursor.fetchall()
text = 'alignmentPosition\thmmPosition\tResidues_at_the_given_position\n'
for align_row in alignments:
    alnpos, pfampos, residues = align_row
    text += str(alnpos) + '\t'
    text += str(pfampos) + '\t'
    text += residues + '\n'

open('alignment_table.tsv', 'w').write(text)
'''
# make ATP table
mycursor.execute("SELECT * FROM ligands where ligand~%s", ('ATP',))
hits = mycursor.fetchall()
text = 'ligand\tligand_id\tWTpos\tAlignmentPos\tUniProtAcc\tGene\tPDB\tPubMed\n'
for hit in hits:
    row_id, ligand, ligand_id, uniprotpos, pfampos, acc, gene, pdb, pubmed = hit
    for item in [ligand, ligand_id, uniprotpos, pfampos, acc, gene, pdb, pubmed]:
        text += str(item).lstrip().rstrip() + '\t'
    text += '\n'

open('atp_table2.tsv', 'w').write(text)