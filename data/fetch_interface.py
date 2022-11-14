#!/usr/bin/env python3
# coding: utf-8

'''
A script to fetch interfaces from the PDB
graph database via API for given accessions
'''

import Bio.PDB, requests
import gzip, os, sys, threading

api_url = 'https://www.ebi.ac.uk/pdbe/graph-api/uniprot/interface_residues/'
LIGANDS = ['ATP', 'ADP', 'MG', 'MN',
            '0WM','1LT','07J','DB8',
            '6GY','4MK','6T2','VGH',
            'P06','1N1','AQ4','E53',
            'IRE','STI','NIL','YY3',
            'LQQ','P30','BAX','B49',
            '032'] ## list of ligand codes to search
PFAM_DOMS = ['PK_Tyr_Ser-Thr', 'Pkinase'] # Domains to search for
PATH_TO_FASTA = '../KA/UniProtFasta2/'
OUT_TEXT = '#acc\tGene\tPfam_Dom\Interface\tPfam_Position\tUniProt_Position\tUniProt_AA\tPDB\n' # Header of the output file
OUT_FILE = 'interface_sites.tsv' # name of output file

class Kinases:
    '''
    Class to define kinases
    '''
    def __init__(self, acc, gene) -> None:
        self.acc = acc
        self.gene = gene
        self.interface = []
        self.kinase_to_pfam = {}
    def append_interface(self, index, code, pdb_entries):
        self.interface.append(BindingSite(index, code, pdb_entries))

class BindingSite:
    '''
    Class to define a ligand binding site
    '''
    def __init__(self, index, code, pdb_entries) -> None:
        self.index = int(index)
        self.code = code
        self.pdb_entries = pdb_entries

'''
read UniProtFasta2 sequences and fetch their ligand
binding information from PDB EBI graph library. It
checks only for those ligands that are specified above
'''
kinase_dic = {}
path_to_fasta = '../KA/UniProtFasta2/'
count = 0
for file in os.listdir(path_to_fasta):
    # count += 1
    # if (count == 50):
    #     break
    # read only FASTA files
    if file.endswith('.fasta') is False:
        continue
    # fetch gene name also from the FASTA file
    for line in open(path_to_fasta+file, 'r'):
        if line[0] == '>':
            gene = line.split('GN=')[1].split()[0]
            break
    # make accession from file name
    acc = file.split('.')[0]
    # Make an instance of the kinase with gene name and acc
    kinase_dic[acc] = Kinases(acc, gene)
    # Fetch response from EBI PDBe graph library
    response = requests.get(api_url+acc)
    if int(response.status_code) != 200:
        # raise ValueError(acc, 'failed')
        print ('Failed to fetch information for', acc, 'from PDBe graph library')
        flag = 0
        for line in open(PATH_TO_FASTA+acc+'.txt', 'r'):
            if line[:2] != 'FT':
                continue
            if line.split()[1] in ['SITE', 'ACT_SITE']:
                # print (line)
                sites = line.split()[2].replace('\n', '')
                if '..' in sites:
                    sites = [int(i) for i in range(int(sites.split('..')[0]), int(sites.split('..')[1])+1)]
                else:
                    sites = [int(sites)]
                flag = 1
            if '/note' in line.split()[1] and flag == 1:
                flag = 0
                note = line.split('/note=')[1].replace('\n', '').replace(' ', '').replace('"', '')
                # print (note, sites)
                for site in sites:
                    kinase_dic[acc].append_interface(site, 'XXX', [])
        # sys.exit()
        continue
    # Save response as dic
    dic = response.json()
    # Read ligands in the dic
    for interface in dic[acc]['data']:
        for residue in interface['residues']:
            start_index = str(residue['startIndex'])
            end_index = residue['endIndex']
            start_code = residue['startCode']
            allPDBEntries = interface['allPDBEntries']
            kinase_dic[acc].append_interface(start_index, start_code, allPDBEntries)
            # print (acc, allPDBEntries)
            # print (acc, ligand['accession'], start_index, start_code, ';'.join(allPDBEntries))
            # OUT_TEXT += acc +'\t'+ ligand['accession'] +'\t'+ start_index +'\t'+ start_code +'\t'+ ';'.join(allPDBEntries) + '\n'
# sys.exit()
# Map ligand binding site back to the
# Pfam domain either Pkinase or Tyr_pkinase
for pfam_dom in PFAM_DOMS:
    HMMSEARCH_OUT = 'allKinasesHmmsearch'+pfam_dom+'.txt'
    flag = 0
    pfam = {}
    ## Read HMMSearch file and save the mapping for a given pfam_domain
    for line in open(HMMSEARCH_OUT, 'r'):
        if len(line.split()) == 0:
            continue
        if line[:2] == '>>':
            flag = 0
            ## lines with kinase start
            kinase = line.split('>>')[1].lstrip().rstrip()
            # Raise an error if the kinase instance not found
            if kinase not in kinase_dic:
                raise ValueError(f'{kinase} not found in the HMMSearch output of {pfam_dom}')
            if pfam_dom not in kinase_dic[kinase].kinase_to_pfam:
                kinase_dic[kinase].kinase_to_pfam[pfam_dom] = {}
            flag = 1
        elif line.split()[0] == pfam_dom:
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
                        kinase_dic[kinase].kinase_to_pfam[pfam_dom][kin_start] = pfam_start
                        pfam_start += 1
                        kin_start += 1
                    elif pfam_char in ['.', '-']:
                        kin_start += 1
                    elif kin_char in ['.', '-']:
                        pfam_start += 1
                    else:
                        print ('Exception found', kinase)
                        sys.exit()
    # print (kinase_dic['P21802'].kinase_to_pfam[628])
    # print (kinase_dic[kinase].kinase_to_pfam)

    # Write the ligand binding information based on the pfam_domain
    for kinase in kinase_dic:
        if pfam_dom not in kinase_dic[kinase].kinase_to_pfam:
            continue
        for binding_site in kinase_dic[kinase].interface:
            # print (binding_site.index)
            if binding_site.index not in kinase_dic[kinase].kinase_to_pfam[pfam_dom]:
                print (kinase, binding_site.index, kinase_dic[kinase].kinase_to_pfam[pfam_dom])
                continue
            OUT_TEXT += kinase +'\t'
            OUT_TEXT += kinase_dic[kinase].gene +'\t'
            OUT_TEXT += pfam_dom +'\t'
            OUT_TEXT += 'INTERFACE' +'\t'
            OUT_TEXT += str(kinase_dic[kinase].kinase_to_pfam[pfam_dom][binding_site.index]) + '\t'
            OUT_TEXT += str(binding_site.index) +'\t'
            OUT_TEXT += binding_site.code +'\t'
            OUT_TEXT += ';'.join(binding_site.pdb_entries) + '\n'

# Save the output file
open(OUT_FILE, 'w').write(OUT_TEXT)