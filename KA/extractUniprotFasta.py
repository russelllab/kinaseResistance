#!/usr/bin/env python3
# coding: utf-8

# A script to convert gene symbols to uniprot accs,
# and extract FASTA formatted sequences

import os,sys
import argparse, requests, urllib

fasta = {}
#for line in open('../fasta_cosmic_kinases_for_alignment_full.fasta', 'r'):
for files in ['kinases_cosmic_uniprot.csv', 'kinases_deact_mutations_uniprot.csv', 'kinases_act_mutations_uniprot.csv']:
    for line in open(files, 'r'):
        #if line[0] == '>':
        if line[0] != '#' and 'uniprot_id' not in line:
            print (files, line)
            #name = line.split('>')[1].split('\n', '')
            #name = 'ALK'
            if 'cosmic' in files:
                acc = line.split(',')[1].replace('\n', '')
                name = line.split(',')[0].replace('\n', '')
            else:
                acc = line.split(',')[0].replace('\n', '')
                name = line.split(',')[1].replace('\n', '')
            if acc not in fasta:
                '''
                url = 'https://www.uniprot.org/uploadlists/'
                #print (name)
                params = {
                'from': 'GENENAME',
                'to': 'ACC',
                'format': 'tab',
                'taxon': '9606',
                'query': name
                }

                data = urllib.parse.urlencode(params)
                data = data.encode('utf-8')
                req = urllib.request.Request(url, data)
                with urllib.request.urlopen(req) as f:
                   response = f.read()
                acc = None
                seq = None
                for line in response.decode('utf-8').split('\n'):
                    if 'From' not in line:
                        if len(line.split())>0:
                            acc = str(line.split('\t')[1].replace('\n', ''))
                            break
                '''
                print (acc)
                #sys.exit()
                if acc != None:
                    response = requests.get('https://www.uniprot.org/uniprot/'+acc+'.fasta')
                    if response.status_code == 200:
                        response = urllib.request.urlopen('https://www.uniprot.org/uniprot/'+acc+'.fasta')
                        seq = '>'+name+'_'+acc+'\n'
                        #seq = ''
                        for line in response.read().decode('utf-8').split('\n'):
                            if len(line.split())>0:
                                if line[0] != '>':
                                    #seq += line.replace('\n', '')
                                    seq += line
                        fasta[acc] = seq
                    else:
                        print ('UniProt Acc not found for', name)

for acc in fasta:
    open('UniProtFasta2/'+acc+'.fasta', 'w').write(fasta[acc])
