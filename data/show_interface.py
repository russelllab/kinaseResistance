#!/usr/bin/python3

'''
This script maps drug resistance mutations in kinases
to 3D structure available from PDB
'''

from Bio.PDB import *
import os, sys, gzip
from pymol import cmd

def plot(acc, pdb, givenChain, inhibitor):
    if os.path.isfile(acc+'.fasta') == False:
        os.system('wget https://www.uniprot.org/uniprot/'+acc+'.fasta')

    l = ''
    for line in open(pdb+'_'+givenChain+'.fasta', 'r'):
        l += line
    l += '\n'
    for line in open(acc+'.fasta', 'r'):
        l += line
    #os.system('cat '+pdb+'_'+givenChain+'.fasta > input.fasta')
    #os.system('cat '+acc+'.fasta >> input.fasta')
    open('input.fasta', 'w').write(l)
    os.system('clustalo -i input.fasta --outfmt=clu --force -o output.txt')

    fastaPDB=[]
    fastaUniprot=[]
    for line in open('output.txt', 'r'):
        if 'CLUSTAL' not in line and line.split()!=[] and line.split()!=['//']:
            #print(line.split())
            if len(line.split()) > 1:
                if pdb in line.split()[0]:
                    fastaPDB += line.replace('\n','').split()[1]
                elif acc in line.split()[0]:
                    fastaUniprot += line.replace('\n','').split()[1]

    #print (fastaUniprot)
    #print (fastaPDB)
    #sys.exit()
    uniprotToFasta = {}
    countPDB=0; countUniprot=0
    for residuePDB, residueUniprot in zip(fastaPDB, fastaUniprot):
        #print (residuePDB, residueUniprot)
        if residuePDB != '-' and residueUniprot != '-':
            uniprotToFasta[str(countUniprot)] = str(countPDB)
            countPDB += 1
            countUniprot += 1
        elif residuePDB == '-':
            countUniprot += 1
        else:
            countPDB += 1

    #print (uniprotToFasta)
    #sys.exit()

    parser = MMCIFParser()
    structure = parser.get_structure(pdb, pdb+".cif")
    fastaToAtom = {}
    for model in structure:
        for chain in model:
            if (chain.id == givenChain):
                fasta = '>'+pdb+'_'+chain.id+'\n'
                map = []; num = 0
                for residue in chain:
                    for atom in residue:
                        if (atom.id == 'CA'):
                            if residue.id[0] == ' ':
                                #print (pdb, chain.id, residue.get_resname(), residue.id)
                                AA = protein_letters_3to1[residue.get_resname()]
                                fasta += AA
                                fastaToAtom[str(num)] = str(residue.id[1])
                                num += 1
                #print (fasta)
                #open(pdb+'_'+chain.id+'.fasta', 'w').write(fasta)
                break

    #print (fastaToAtom)
    #sys.exit()

    interfaces = []
    for line in open('interface.tsv', 'r'):
        if line.split('\t')[0] == pdb and line.split('\t')[1] == givenChain and line.split('\t')[4] == inhibitor:
            interfaces.append(str(line.split('\t')[3]))

    resitancePositions = []
    for line in open('../KA/FINAL_FILE_MUTATIONS.txt', 'r'):
        if 'Original.mech' not in line:
            mechismoFormat = line.split(',')[2].replace('"', '')
            uniprot = mechismoFormat.split('/')[0]
            mutation = mechismoFormat.split('/')[1]
            mutationPosition = str(int(mutation[1:-1]) - 1)
            if uniprot == acc:
                if mutationPosition in uniprotToFasta:
                    fastaPosition = uniprotToFasta[mutationPosition]
                    if fastaPosition in fastaToAtom:
                        atomPosition = fastaToAtom[fastaPosition]
                        resitancePositions.append(str(atomPosition))
                        print (mutationPosition, fastaPosition, atomPosition)

    #sys.exit()
    print (resitancePositions)
    print (interfaces)
    print (fastaToAtom)
    print (uniprotToFasta)
    #sys.exit()

    cmd.load(pdb+'.cif')
    cmd.color('grey', 'chain '+givenChain)
    cmd.color('yellow', 'resn '+inhibitor)
    cmd.hide('everything')
    cmd.show('cartoon', 'chain '+str(givenChain))
    cmd.show('licorice', 'resn '+str(inhibitor))
    #cmd.remove('resn hoh')
    cmd.set('cartoon_fancy_helices', 1)
    cmd.bg_color('white')
    cmd.set ("sphere_scale", 0.5)

    for position in resitancePositions:
        if position not in interfaces:
            cmd.color('red', givenChain+'/'+position+'/ca')
            cmd.show('spheres', givenChain+'/'+position+'/ca')
        else:
            cmd.color('cyan', givenChain+'/'+position+'/ca')
            cmd.show('spheres', givenChain+'/'+position+'/ca')

    for position in interfaces:
        if position not in resitancePositions:
            cmd.color('green', givenChain+'/'+position+'/ca')
            cmd.show('spheres', givenChain+'/'+position+'/ca')

    '''
    for line in open('interface.tsv', 'r'):
        if line.split('\t')[0] == pdb and line.split('\t')[1] == givenChain and line.split('\t')[4] == inhibitor:
            position = str(line.split('\t')[3])
            cmd.color('green', givenChain+'/'+position+'/ca')
            cmd.show('spheres', givenChain+'/'+position+'/ca')
    '''

pdb = '3ug2'
givenChain = 'A'
inhibitor = 'IRE'
acc = 'P00533'
plot(acc, pdb, givenChain, inhibitor)
