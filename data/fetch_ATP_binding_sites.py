#!/usr/bin/env python3
# coding: utf-8

'''
A script to fetch ATP binding sites from
PDB graph database via API for given accessions
'''

import Bio.PDB, requests
import gzip, os, sys, threading

api_url = 'https://www.ebi.ac.uk/pdbe/graph-api/uniprot/ligand_sites/'
LIGANDS = ['ATP', 'ADP', 'MG', 'MN'] ## list of ligand codes to search
PFAM_DOMS = ['PK_Tyr_Ser-Thr', 'Pkinase'] # Domains to search for
OUT_TEXT = '#acc\tGene\tPfam_Dom\tLigand\tPfam_Position\tUniProt_Position\tUniProt_AA\tPDB\n' # Header of the output file
OUT_FILE = 'ATP_binding_sites3.tsv' # name of output file

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

'''
read UniProtFasta2 sequences and fetch their ligand
binding information from PDB EBI graph library. It
checks only for those ligands that are specified above
'''
kinase_dic = {}
path_to_fasta = '../KA/UniProtFasta2/'
for file in os.listdir(path_to_fasta):
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
        continue
    # Save response as dic
    dic = response.json()
    # Read ligands in the dic
    for ligand in dic[acc]['data']:
        # consider only the specified ligands
        if ligand['accession'] not in LIGANDS:
            continue
        for residue in ligand['residues']:
            start_index = str(residue['startIndex'])
            end_index = residue['endIndex']
            start_code = residue['startCode']
            allPDBEntries = residue['allPDBEntries']
            kinase_dic[acc].append_ligand(ligand['accession'], start_index, start_code, allPDBEntries)
            # print (acc, ligand['accession'], start_index, start_code, ';'.join(allPDBEntries))
            # OUT_TEXT += acc +'\t'+ ligand['accession'] +'\t'+ start_index +'\t'+ start_code +'\t'+ ';'.join(allPDBEntries) + '\n'

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
            ## lines with kinase start
            kinase = line.split('>>')[1].lstrip().rstrip()
            # Raise an error if the kinase instance not found
            if kinase not in kinase_dic:
                raise ValueError(f'{kinase} not found in the HMMSearch output')
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
    print (kinase_dic[kinase].kinase_to_pfam)

    # Write the ligand binding information based on the pfam_domain
    for kinase in kinase_dic:
        if pfam_dom not in kinase_dic[kinase].kinase_to_pfam:
            continue
        for binding_site in kinase_dic[kinase].ligand_binding_sites:
            # print (binding_site.index)
            if binding_site.index not in kinase_dic[kinase].kinase_to_pfam[pfam_dom]:
                print (kinase, binding_site.index, kinase_dic[kinase].kinase_to_pfam[pfam_dom])
                continue
            OUT_TEXT += kinase +'\t'
            OUT_TEXT += kinase_dic[kinase].gene +'\t'
            OUT_TEXT += pfam_dom +'\t'
            OUT_TEXT += binding_site.ligand +'\t'
            OUT_TEXT += str(kinase_dic[kinase].kinase_to_pfam[pfam_dom][binding_site.index]) + '\t'
            OUT_TEXT += str(binding_site.index) +'\t'
            OUT_TEXT += binding_site.code +'\t'
            OUT_TEXT += ';'.join(binding_site.pdb_entries) + '\n'

# Save the output file
open(OUT_FILE, 'w').write(OUT_TEXT)