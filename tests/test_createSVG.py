#!/usr/bin/env python3.10
# -*- coding: utf-8 -*-
'''
A pytest script to test create_svg_2023*_kinases*.py
'''

import pytest
import psycopg2
import os, sys, ast
os.sys.path.append('ML/')
import fetchData
sys.path.insert(1, 'Create_SVG/Enhancements_May2023/')
import create_svg_20230510_kinases_GS as create_svg
conservation_dic_path = 'Create_SVG/Enhancements_May2023/'+'GenerelleKonservierung_May-10-2023.txt'
identity_dic_path = 'Create_SVG/Enhancements_May2023/'+'SeqIdentity_Matrix_May-10-2023.txt'

def get_accs_in_alignment(alignmentFile):
    """
    Get the accessions in the alignment file
    """
    accs_in_alignment = {}
    with open(alignmentFile) as f:
        for line in f:
            acc = line.split('|')[1]
            name = line.split()[0]
            accs_in_alignment[acc] = name
    return accs_in_alignment

def make_dic_mutation_info(accs_in_alignment):
    """
    Make a dictionary of kinase mutations in the alignment
    by connecting to the DB and querying the mutations table
    """
    # Connect to the DB
    mydb = fetchData.connection()
    mydb.autocommit = True
    mycursor = mydb.cursor()
    # Get the mutations
    mycursor.execute("SELECT wtpos, mut_type, acc FROM mutations")
    hits = mycursor.fetchall()
    # Make the dictionary of to map mutations types
    # to their full names
    mut_types = {'A': 'Activating', 'D': 'Deactivating', 'R': 'Resistance'}
    # Make the dictionary of mutations info
    dic_mutations_info = {}
    # Query the hits
    for hit in hits:
        position, mutType, acc = hit
        # Ignore the acc not in the alignment
        if acc not in accs_in_alignment: continue
        # Ignore the neutral mutations
        if mutType == 'N': continue
        if acc not in dic_mutations_info:
            dic_mutations_info[acc] = {'Activating':[], 'Deactivating':[], 'Resistance':[]}
        # extract name of the sequence
        # given the accession
        name = accs_in_alignment[acc]
        # extract the start position
        start = int(name.split('|')[-1])
        # Ignore the mutations before the start position
        if int(position) >= start:
            dic_mutations_info[acc][mut_types[mutType]].append(str(position))
    return dic_mutations_info
    
def test_createSVG():
    """
    Test the CreateSVG script
    """
    with open(conservation_dic_path) as g:
        overconserv = g.read()
        overallconservation = ast.literal_eval(overconserv)
    with open(identity_dic_path) as h:
        data_ident = h.read()
    identitydictionary = ast.literal_eval(data_ident)
    alignmentFile = 'webApp/static/hmm/humanKinasesTrimmed.clustal'
    accs_in_alignment = get_accs_in_alignment(alignmentFile)
    # Instances to test
    instances = [
                    ['P15056', '600', True],
                    ['P1506', '600', False], # wrong acc, should throw a KeyError
                    ['O96017', '373', True],
                ]
    for instance in instances:
        for ws in [10, 25, 50]:
            for topN in [10, 20, 50]:
                for sortingvalue in ['1', '2']:
                    # print (instance[0], ws, topN)
                    known_verdict = instance[2]
                    mutation_position = int(instance[1])
                    kinase = instance[0]
                    dic_mutations_info = make_dic_mutation_info(accs_in_alignment)
                    verdict = kinase in dic_mutations_info
                    # check if the verdict is correct (as specified)
                    assert verdict == known_verdict
                    # run the create_svg script
                    try:
                        create_svg.main(sortingvalue, identitydictionary, overallconservation,\
                                        alignmentFile, kinase, mutation_position,\
                                        int(ws), int(topN), dic_mutations_info)
                    except KeyError as error:
                        # assert if this happens when the given kinase is in the list below
                        assert kinase in ['P1506']
                        continue

if __name__ == '__main__':
    test_createSVG()
