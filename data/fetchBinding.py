#!/usr/bin/env python3
# coding: utf-8

'''
read .txt file from UniProt for a protein and extract the
BINDING section information and save it in a dictionary
'''

import os, sys, gzip
import re
import pandas as pd
import argparse

class Binding:
    '''
    class to store binding information
    '''
    def __init__(self, position):
        self.position = position
        self.text = ''

    def extractInformation(self):
        '''
        extract ligand, ligand_id, evidence and pubmed information from the text
        '''
        text = self.text
        note_pattern = r'/note="([^"]*)'
        ligand_pattern = r'/ligand="([^"]*)'
        ligand_id_pattern = r'/ligand_id="([^"]*)"'
        evidence_pattern = r'/evidence="([^"]*)"'
        pdb_pattern = r"PDB:([A-Z0-9]+)"
        pubmed_pattern = r"PubMed:\s*(\d+)"

        note_matches = re.findall(note_pattern, text)
        ligand_matches = re.findall(ligand_pattern, text)
        ligand_id_matches = re.findall(ligand_id_pattern, text)
        evidence_matches = re.findall(evidence_pattern, text)
        pdb_matches = re.findall(pdb_pattern, text)
        pubmed_matches = re.findall(pubmed_pattern, text)

        note = ''.join([match for match in note_matches])
        ligand = ','.join([match for match in ligand_matches])
        ligand_id = ','.join([match for match in ligand_id_matches])
        evidence = ''.join([match for match in evidence_matches])
        pdb = ','.join([match for match in pdb_matches])
        pubmed = ','.join([match for match in pubmed_matches])

        return ligand, ligand_id, evidence, pdb, pubmed

def is_uppercase_with_special_chars(string):
    if string.isupper() and re.match(r'^[A-Z_]+$', string):
        return True
    else:
        return False

def main(file):
    '''
    read .txt file from UniProt for a protein and extract the
    BINDING section information and save it in a dataframe.
    '''
    acc, gene = None, None
    flag = 0
    dic_binding = {}
    # for line in gzip.open('../KA/UniProtFasta2/'+acc+'.txt.gz', 'rt'):
    for line in gzip.open(file, 'rt'):
        if line.startswith('AC'):
            if acc is None: acc = line.split()[1].split()[0].replace(';', '')
            continue
        if line.startswith('GN'):
            if gene is None: gene = line.split('Name=')[1].split()[0].replace(';', '')
            continue
        # print (line[3:12])
        if line.startswith('FT') is False: continue
        if line.startswith('FT   BINDING'):
            # check if the binding site is number
            sites = line.split()[2].rstrip()
            if '..' in sites:
                start, end = sites.split('..')
                if start.isdigit() is False or end.isdigit() is False:
                    flag = 0
                    continue
                start, end = int(start), int(end)
                sites = range(start, end+1)
            elif sites.isdigit() is False:
                print ('Error:', acc, gene, line)
                flag = 0
                continue
            else:
                sites = [int(sites)]
            # print (position)
            flag = 1
            continue
        # elif line.split()[1].isupper():
        #     flag = 0
        else:
            if is_uppercase_with_special_chars(line.split()[1]) is True:
                flag = 0
                continue

        if flag != 1: continue
        for site in sites:
            if site not in dic_binding:
                dic_binding[site] = Binding(site)
            dic_binding[site].text += line.replace('\n', ' ').lstrip('FT').lstrip(' ')

    data = []
    for site in dic_binding:
        # print (mutation, dic_mutations[mutation].text)
        ligand, ligand_id, evidence, pdb, pubmed = dic_binding[site].extractInformation()
        if ligand_id == '': ligand_id = '-'
        if evidence == '': evidence = '-'
        if pdb == '': pdb = '-'
        if pubmed == '': pubmed = '-'
        # print (acc, gene, mutation, dic_mutations[mutation].mut_type, note, evidence, pubmed)
        row= []
        row.append(acc)
        row.append(gene)
        row.append(site)
        row.append(ligand)
        row.append(ligand_id)
        row.append(evidence)
        row.append(pdb)
        row.append(pubmed)
        data.append(row)

    df = pd.DataFrame(data, columns = ['UniProt',
                                    'Gene',
                                    'Site',
                                    'Ligand',
                                    'LigandID',
                                    'Evidence',
                                    'PDB',
                                    'Pubmed']
                                    )
    print (df.to_csv(index=False, sep='\t'))
    return (df.to_csv(index=False, sep='\t'))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Fetch variants/mutagens from UniProt file')
    parser.add_argument('file', help='path to UniProt file')
    args = parser.parse_args()
    file = args.file
    main(file)
    # acc = 'P00533'
    # acc = 'Q13546'
    # acc = 'P22607'