#!/usr/bin/env python3
# coding: utf-8

import Bio.PDB
import gzip

pdbs = []
for line in gzip.open('PF07714_pdb_chain_pfam.tsv.gz', 'rt'):
    pdbs.append(line.split('\t')[0])

pdbs = list(set(pdbs))

for pdb in pdbs:
    dic[pdb] = Bio.PDB.MMCIF2Dict.MMCIF2Dict('3ug2.cif')
    print(dic[pdb]['_pdbx_entity_nonpoly.name'])
    print (dic[pdb]['_pdbx_entity_nonpoly.comp_id'])
