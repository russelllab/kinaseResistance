from Bio.Data.IUPACData import protein_letters_3to1
from Bio.PDB import *
import numpy as np
import os, sys

pdbs = {}
## Load inhibitors from HETATM file ############################################
class inhibitors:
    def __init__(self, code, name, pdb):
        self.name = name
        self.code = code
        self.pdb = [pdb]

inhibitorsDic = {}
for line in open('HETATM.tsv', 'r'):
    if line[0] != '#':
        inhibitorName = line.replace('\n', '').split('\t')[2]
        inhibitorCode = line.split('\t')[1]
        inhibitorPDB = line.split('\t')[0]
        pdbs[inhibitorPDB] = []
        if 'inib' in inhibitorName:
            if inhibitorCode not in inhibitorsDic:
                inhibitorsDic[inhibitorCode] = inhibitors(inhibitorCode, inhibitorName, inhibitorPDB)
            else:
                inhibitorsDic[inhibitorCode].pdb.append(inhibitorPDB)
########################################################################################


## Load chain with kinases from filtered SIFTS file
for line in open('PF07714_pdb_chain_pfam.tsv', 'r'):
    pdb = line.split('\t')[0]
    chain = line.split('\t')[1]
    if pdb in pdbs:
        pdbs[pdb].append(chain)
########################################################################################


def run(pdb, given_chain, given_ligand, ligand_name, l):
    tempdir = ''
    #pdb = '5MO4'
    format = 'cif'
    #given_chain = 'A'
    #given_ligand = 'NIL'

    ## 3 letters to 1 letter dictionary of Amino Acids
    AA ={'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q',
    'ASP':'D', 'ASN':'N', 'HIS':'H', 'TRP':'W', 'PHE':'F', 'TYR':'Y',
    'ARG':'R', 'LYS':'K', 'SER':'S', 'THR':'T', 'MET':'M', 'ALA':'A',
    'GLY':'G', 'PRO':'P', 'CYS':'C'}

    ## To check if both chains are present in the complex
    if format == 'cif':
        try:
            parser = MMCIFParser()
            structure = parser.get_structure(pdb, tempdir+pdb+'.cif')
        except:
            os.system('wget https://files.rcsb.org/view/'+pdb+'.cif -nv -O '+tempdir+pdb+'.cif')
            parser = MMCIFParser()
            structure = parser.get_structure(pdb, tempdir+pdb+'.cif')
    else:
        try:
            parser = PDBParser()
            structure = parser.get_structure(pdb, tempdir+pdb+'.pdb')
        except:
            os.system('wget https://files.rcsb.org/view/'+pdb+'.pdb -nv -O '+tempdir+pdb+'.pdb')
            parser = PDBParser()
            structure = parser.get_structure(pdb, tempdir+pdb+'.pdb')

    for model in structure:
        if given_chain in model:
            chain = model[given_chain]
        else:
            print (f'Chain {given_chain} absent in {pdb}')
            sys.exit()

        fasta1 = ''
        for count, residue1 in enumerate(chain):
            het1 = residue1.get_id()[0]
            #print (len(residue1.get_id()[0]))
            if het1 == ' ':
                #print (residue1.get_resname())
                if residue1.get_resname() in AA:
                    fasta1 += protein_letters_3to1[residue1.get_resname()]
                    for atom1 in residue1:
                        for residue2 in chain:
                            het2 = residue2.get_id()[0]
                            if het2 != ' ':
                                #print (residue2.get_resname())
                                if residue2.get_resname() == given_ligand:
                                    for atom2 in residue2:
                                        diff_VdW  = atom1.coord - atom2.coord
                                        if np.sqrt(np.sum(diff_VdW * diff_VdW)) <= 6.5:
                                            distance = np.sqrt(np.sum(diff_VdW * diff_VdW))
                                            #print (chain.id, residue1.get_resname(), residue1.id[1], residue2.get_resname(), residue2.id[1], distance)

                                            l += pdb +'\t'+\
                                            chain.id +'\t'+\
                                            protein_letters_3to1[residue1.get_resname()] +'\t'+\
                                            str(residue1.id[1]) +'\t'+\
                                            residue2.get_resname() +'\t'+\
                                            ligand_name +'\t'+\
                                            str(residue2.id[1]) +'\t'+\
                                            str(distance) + '\n'
                                            #closest_residue = residueA2.id[1]

        print (fasta1)
        open(pdb+'_'+given_chain+'.fasta', 'w').write('>'+pdb+'|'+given_chain+'\n'+fasta1)
        #sys.exit()
    return l

l = '#PDB\tChain\tAA\taaKin\tposKin\tcodeInhib\tnameInhib\tposInhib\tDistance\n'
for inhibitorCode in inhibitorsDic:
    for pdb in inhibitorsDic[inhibitorCode].pdb:
        for chain in pdbs[pdb]:
            l = run(pdb, chain, inhibitorCode, inhibitorsDic[inhibitorCode].name, l)
            #sys.exit()
open('interface.tsv', 'w').write(l)
