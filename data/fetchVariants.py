#!/usr/bin/env python3
# coding: utf-8

'''
read .txt file from UniProt for a protein and extract the
VARIANT/MUTAGENESIS section information and save it in a dictionary
'''

import os, sys, gzip
import re
import pandas as pd
import argparse

class Mutation:
    '''
    class to store mutation information
    '''
    def __init__(self, mutation, mut_type):
        self.mutation = mutation
        self.mut_type = mut_type
        self.text = ''

    def extractInformation(self):
        '''
        extract note, evidence and pubmed information from the text
        '''
        text = self.text
        note_pattern = r'/note="([^"]*)'
        evidence_pattern = r'/evidence="([^"]*)"'
        pubmed_pattern = r"PubMed:\s*(\d+)"

        note_matches = re.findall(note_pattern, text)
        evidence_matches = re.findall(evidence_pattern, text)
        pubmed_matches = re.findall(pubmed_pattern, text)

        note = ''.join([match for match in note_matches])
        evidence = ''.join([match for match in evidence_matches])
        pubmed = ','.join([match for match in pubmed_matches])

        return note, evidence, pubmed

def main(file):
    '''
    read .txt file from UniProt for a protein and extract the
    VARIANT/MUTAGEN section information and save it in a 
    dataframe.
    '''
    acc, gene = None, None
    flag = 0
    dic_mutations = {}
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
        if line.startswith('FT   VARIANT') or line.startswith('FT   MUTAGEN'):
            # check if the variant psoition is number
            position = line.split()[2].rstrip()
            if not position.isdigit():
                flag = 0
                continue
            # print (position)
            flag = 1
            note = 0
            mut_type = line.split()[1]
            continue    
        elif line.split()[1].isupper() and line.split()[1].isalpha():
            flag = 0

        if flag != 1: continue
        if '/note=' in line:
            #pattern to extract mutation
            mutation_pattern = r'("\S+)\s*->\s*(\S+)'
            # mutation_pattern = r'note="([\w\s]+)->([\w\s,]+?)(?:, | or )?([\w\s]+)?"'
            mutation_text = re.search(mutation_pattern, line)
            try:
                mutation_text = mutation_text.group().replace('"', '').replace(':', '').replace(' ', '')
            except AttributeError as e:
                print ('Error:', acc, gene, e, line)
                flag = 0
                note = 0
                continue
            wild_type = mutation_text.split('->')[0]
            if ',' in mutation_text:
                mutants = mutation_text.split('->')[1].split(',')
            elif 'or' in mutation_text:
                mutants = mutation_text.split('->')[1].split('or')
            else:
                if mutation_text[0].isalpha() and mutation_text[-1].isalpha() and len(mutation_text.replace('->','')) == 2:
                    mutants = [mutation_text[-1]] 
                else:
                    mutants = [mutation_text.split('->')[1]]
                    # print ('Error:', acc, gene, line)
                    # flag = 0
                    # note = 0
                    # continue
            for mutant in mutants:
                mutation = wild_type + position + mutant
                mutation = mutation.lstrip().rstrip()
                if mutation not in dic_mutations:
                    dic_mutations[mutation] = Mutation(mutation, mut_type)
                dic_mutations[mutation].text = line.replace('\n', ' ').lstrip('FT').lstrip(' ')
                # print (position, mutant, mutation, dic_mutations[mutation].text)
            note = 1
            # print (mutation, dic_mutations[mutation].text)
        elif note == 1:
            # dic_mutations[mutation].text += line.lstrip('FT').lstrip(' ').replace('\n', ' ')
            for mutant in mutants:
                mutation = wild_type + position + mutant
                mutation = mutation.lstrip().rstrip()
                dic_mutations[mutation].text += line.replace('\n', ' ').lstrip('FT').lstrip(' ')

    data = []
    for mutation in dic_mutations:
        # print (mutation, dic_mutations[mutation].text)
        note, evidence, pubmed = dic_mutations[mutation].extractInformation()
        # print (acc, gene, mutation, dic_mutations[mutation].mut_type, note, evidence, pubmed)
        row= []
        row.append(acc)
        row.append(gene)
        row.append(mutation)
        row.append (acc + '/' + mutation)
        row.append (gene + '/' + mutation)
        row.append(dic_mutations[mutation].mut_type)
        row.append(note)
        row.append(evidence)
        row.append(pubmed)
        data.append(row)

    df = pd.DataFrame(data, columns = ['UniProt',
                                    'Gene',
                                    'Mutation',
                                        'UniProt/Mutation',
                                        'Gene/Mutation',
                                    'Type',
                                    'Note',
                                    'Evidence',
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