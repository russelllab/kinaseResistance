#!/usr/bin/env python3
# coding: utf-8

'''
A script to fetch ATP binding sites from
PDB graph database via API for given accessions
'''

import Bio.PDB, requests
import gzip, os, sys, threading

api_url = 'https://www.ebi.ac.uk/pdbe/graph-api/uniprot/ligand_sites/'

class Kinases:
    '''
    Class to define kinases
    '''
    def __init__(self, acc, gene) -> None:
        self.acc = acc
        self.gene = gene
        self.ligand_binding_sites = []
        self.kinase_to_pfam = {}
    def append_ligand(self, ligand, index, code, pdb_entries):
        self.ligand_binding_sites.append(BindingSite(index, code, ligand, pdb_entries))

class BindingSite:
    '''
    Class to define a ligand binding site
    '''
    def __init__(self, index, code, ligand, pdb_entries) -> None:
        self.index = int(index)
        self.code = code
        self.ligand = ligand
        self.pdb_entries = pdb_entries

## Store accessions
kinase_dic = {}
path_to_fasta = '../KA/UniProtFasta2/'
for file in os.listdir(path_to_fasta):
    if file.endswith('.fasta') is False:
        continue
    for line in open(path_to_fasta+file, 'r'):
        if line[0] == '>':
            gene = line.split('GN=')[1].split()[0]
            break
    acc = file.split('.')[0]
    kinase_dic[acc] = Kinases(acc, gene)
    response = requests.get(api_url+acc)
    if int(response.status_code) != 200:
        # raise ValueError(acc, 'failed')
        print (acc, 'failed')
        continue
    dic = response.json()
    for ligand in dic[acc]['data']:
        if ligand['accession'] not in ['ATP', 'ADP']:
            continue
        for residue in ligand['residues']:
            start_index = str(residue['startIndex'])
            end_index = residue['endIndex']
            start_code = residue['startCode']
            allPDBEntries = residue['allPDBEntries']
            kinase_dic[acc].append_ligand(ligand['accession'], start_index, start_code, allPDBEntries)
            # print (acc, ligand['accession'], start_index, start_code, ';'.join(allPDBEntries))
            # out_text += acc +'\t'+ ligand['accession'] +'\t'+ start_index +'\t'+ start_code +'\t'+ ';'.join(allPDBEntries) + '\n'

## Read HMMSearch file
flag = 0
pfam = {}
for line in open('allKinasesHmmsearch.txt', 'r'):
    if len(line.split()) == 0:
        continue
    if line[:2] == '>>':
        ## lines with kinase start
        kinase = line.split('>>')[1].lstrip().rstrip()
        # Raise an error if the kinase instance not found
        if kinase not in kinase_dic:
            raise ValueError(f'{kinase} not found in the HMMSearch output')
        flag = 1
    elif line.split()[0] == 'Pkinase':
        ## lines with Pkinase domain
        pfam_start, pfam_seq, pfam_end = int(line.split()[1]), line.split()[2], int(line.split()[3])
        count = int(line.split()[1])
        for char in pfam_seq:
            if char not in ['.', '-']:
                pfam[count] = char+str(count)
                count += 1
    elif flag == 1:
        if kinase == line.split()[0]:
            ## lines with kinase
            kin_start, kin_seq, kin_end = int(line.split()[1]), line.split()[2], int(line.split()[3])
            for pfam_char, kin_char in zip(pfam_seq, kin_seq):
                if pfam_char not in ['.', '-'] and kin_char not in ['.', '-']:
                    kinase_dic[kinase].kinase_to_pfam[kin_start] = pfam_start
                    pfam_start += 1
                    kin_start += 1
                elif pfam_char in ['.', '-']:
                    kin_start += 1
                elif kin_char in ['.', '-']:
                    pfam_start += 1
                else:
                    print ('Exception found', kinase)
                    sys.exit()
print (kinase_dic['P21802'].kinase_to_pfam[628])

out_text = '#acc\tGene\tLigand\tPfam_Position\tUniProt_Position\tUniProt_AA\tPDB\n'
for kinase in kinase_dic:
    for binding_site in kinase_dic[kinase].ligand_binding_sites:
        # print (binding_site.index)
        if binding_site.index not in kinase_dic[kinase].kinase_to_pfam:
            print (kinase, binding_site.index, kinase_dic[kinase].kinase_to_pfam)
            continue
        out_text += kinase +'\t'
        out_text += kinase_dic[kinase].gene +'\t'
        out_text += binding_site.ligand +'\t'
        out_text += str(kinase_dic[kinase].kinase_to_pfam[binding_site.index]) + '\t'
        out_text += str(binding_site.index) +'\t'
        out_text += binding_site.code +'\t'
        out_text += ';'.join(binding_site.pdb_entries) + '\n'

open('ATP_binding_sites.tsv', 'w').write(out_text)