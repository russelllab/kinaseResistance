#/usr/bin/env python3
# coding: utf-8

'''
A script to take the SIFT pdb_chain to pfam file
as input and extract all the compound information
for the given PDBs and write it in a file
'''

import Bio.PDB
import gzip, os, sys, threading

## Input
pdb_chain_to_pfam_domain = 'PF07714_pdb_chain_pfam.tsv'
pathToPDB = '/net/home.isilon/ds-russell/pdb-cif/' ## Path with CIF formatted files
outputFile = 'HETATM.tsv' ## Output file
threads = 20 ## Number of threads
verbose = 1 ## Verbosity (Binary)

## Multithreading
class myThread (threading.Thread):
    '''
    Class to perform multi-threading
    '''
    def __init__(self, threadID, pdb):
        threading.Thread.__init__(self)
        self.threadID = threadID
        self.pdb = pdb

    def start(self):
        if verbose == 1:
            print ("Threading", self.threadID)
        parseCIF(self.pdb)
        if verbose == 1:
            print ("Exiting", self.threadID)

def call_main(thread, count, pdb):
    thread[count] = myThread(count, pdb)
    thread[count].start()
    while threading.activeCount() == threads:
        pass

'''
Store all PDB accessions in a list
'''
pdbs = []
for line in open(pdb_chain_to_pfam_domain, 'r'):
    pdbs.append(line.split('\t')[0])
pdbs = list(set(pdbs))
print (len(pdbs))

dic ={}
def parseCIF(pdb):
    '''
    Function to parse the PDB file, and store
    pdb id as well as compound id/name in a dic
    '''
    parser = Bio.PDB.MMCIF2Dict.MMCIF2Dict(pdb+'.cif')
    os.system('rm -rf '+pdb+'.cif')
    if '_pdbx_entity_nonpoly.comp_id' in parser:
        dic[pdb] = {}
        #print(dic[pdb]['_pdbx_entity_nonpoly.name'])
        #print (dic[pdb]['_pdbx_entity_nonpoly.comp_id'])
        for compid, compname in zip(parser['_pdbx_entity_nonpoly.comp_id'], parser['_pdbx_entity_nonpoly.name']):
            #l += pdb + '\t' + compid + '\t' + compname + '\n'
            dic[pdb][compid] = compname

'''
A snippet to create an instance of threading
per PDB id. Each instance calls the function
parseCIF, which stores the compound name/id
in a dic
'''
count = 0
thread = {}
for pdb in pdbs:
    code = pdb[1:-1]
    if os.path.isfile(pathToPDB+code+'/'+pdb+'.cif.gz') == True:
        os.system('cp '+pathToPDB+code+'/'+pdb+'.cif.gz .')
        os.system('gunzip '+pdb+'.cif.gz')
        call_main(thread, count, pdb)
        count += 1
        #if (count == 100):
        #    break

## Wait until all threads are not over
while threading.activeCount() > 1:
    pass

## Save the dictionary information
l = '#PDB\tCompID\tCompName\n'
for pdb in dic:
	for compid in dic[pdb]:
		l += pdb + '\t' + compid + '\t' + dic[pdb][compid] + '\n'
open(outputFile, 'w').write(l)
