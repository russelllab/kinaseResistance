#/usr/bin/env python3
# coding: utf-8

import Bio.PDB
import gzip, os, sys, threading

threads = 20
verbose = 1

## Multithreading
class myThread (threading.Thread):
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

pdbs = []
for line in open('PF07714_pdb_chain_pfam.tsv', 'r'):
    pdbs.append(line.split('\t')[0])

pdbs = list(set(pdbs))
print (len(pdbs))

dic ={}
def parseCIF(pdb):
	parser = Bio.PDB.MMCIF2Dict.MMCIF2Dict(pdb+'.cif')
	os.system('rm -rf '+pdb+'.cif')
	if '_pdbx_entity_nonpoly.comp_id' in parser:
		dic[pdb] = {}
		#print(dic[pdb]['_pdbx_entity_nonpoly.name'])
		#print (dic[pdb]['_pdbx_entity_nonpoly.comp_id'])
		for compid, compname in zip(parser['_pdbx_entity_nonpoly.comp_id'], parser['_pdbx_entity_nonpoly.name']):
			#l += pdb + '\t' + compid + '\t' + compname + '\n'
			dic[pdb][compid] = compname

count = 0
thread = {}
pathToPDB = '/net/home.isilon/ds-russell/pdb-cif/'
for pdb in pdbs:
	code = pdb[1:-1]
	if os.path.isfile(pathToPDB+code+'/'+pdb+'.cif.gz') == True:
                os.system('cp '+pathToPDB+code+'/'+pdb+'.cif.gz .')
                os.system('gunzip '+pdb+'.cif.gz')
                call_main(thread, count, pdb)
                count += 1
                #if (count == 100):
                #    break

while threading.activeCount() > 1:
    pass

l = 'PDB\tCompID\tCompName\n'
for pdb in dic:
	for compid in dic[pdb]:
		l += pdb + '\t' + compid + '\t' + dic[pdb][compid] + '\n'
open('HETATM.tsv', 'w').write(l)
