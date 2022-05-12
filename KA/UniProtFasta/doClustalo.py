#!/usr/bin/env python3
# coding: utf-8

# A script to align AND map COSMIC FASTA sequence
# to UniProt FASTA sequences

import os

class kinases:
    def __init__(self, name):
        self.name = name
        self.cosmicSeq = ''
        self.uniprotSeq = ''
        self.cosmicAln = ''
        self.uniprotAln = ''

kins = {}
## Load COSMIC sequences
flag = 1
for line in open('../fasta_cosmic_kinases_for_alignment_full.fasta', 'r'):
    if line[0] == '>':
        name = line.split('>')[1].replace('\n', '')
        flag = 1
        if '_' in name or '.' in name:
            flag = 0
            continue
        print (name)
        kins[name] = kinases(name)

    elif flag == 1:
        kins[name].cosmicSeq += line.replace('\n', '')

## Load UniProt sequences
for name in kins:
    for line in open(name+'_UniProt.fasta', 'r'):
        if line[0] != '>':
            kins[name].uniprotSeq += line.replace('\n', '')

## Write COSMIC and UniProt sequences, and run clustalo
for name in kins:
    l = ''
    l += '>'+name+'_COSMIC\n'
    l += kins[name].cosmicSeq
    l += '\n>'+name+'_UniProt\n'
    l += kins[name].uniprotSeq
    open(name+'_clustalo_input.fasta', 'w').write(l)
    os.system('clustalo --force -i '+name+'_clustalo_input.fasta -o '+name+'_clustalo.aln')

## Map COSMIC to UniProt sequences
for name in kins:
    for line in open(name+'_clustalo.aln', 'r'):
        if line[0] == '>':
            kinSource = line.split('>')[1].replace('\n', '').split('_')[1]
        else:
            if kinSource == 'COSMIC':
                kins[name].cosmicAln += line.replace('\n', '')
            else:
                kins[name].uniprotAln += line.replace('\n', '')

    map = {}
    numCosmic = 0
    numUniprot = 0
    l = '#CosmicPosition\tCosmicAA\tUniProtPosition\tUniProtAA\n'
    for charCosmic, charUniprot in zip(kins[name].cosmicAln, kins[name].uniprotAln):
        if charCosmic!='-' and charUniprot!='-':
            map[numCosmic] = numUniprot
            l += str(numCosmic) + '\t' + str(charCosmic) + '\t' + str(numUniprot) + '\t' + str(charUniprot) + '\n'
            numCosmic += 1
            numUniprot += 1
        elif charCosmic=='-':
            numUniprot += 1
        else:
            numCosmic += 1

    open(name+'_mappings.tsv', 'w').write(l)
