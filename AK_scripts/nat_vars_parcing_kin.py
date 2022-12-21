#parcing natural variants
import os
import pandas as pd
import shutil
#select from robScores2 the kinase uniprot IDs and save them as a list for parcing in two formats once as
#first 4 letters and another as full
uniprot_kinases_raw = os.listdir('/net/home.isilon/ds-russell/kinaseResistance/KA/robScores2/')
#print(uniprot_kinases_raw)
uniprot_kinases = []
uniprot_kinases_short = []
for uniprot_kinase in uniprot_kinases_raw:
    if uniprot_kinase.endswith('.tsv.gz'):
        uniprot_kinase_id = uniprot_kinase[:-7]
        uniprot_kinases.append(uniprot_kinase_id)
        uniprot_kinase_id_short = uniprot_kinase[:4]
        uniprot_kinases_short.append(uniprot_kinase_id_short)
 
#parce through the folders and extract files with kinase uniprot ids and copy them into another folder
for root, dirs, files in os.walk('/net/home.isilon/ds-russell/mechismoX/analysis/mutations/data/VLatest/'):
    dirs[:] = [d for d in dirs if d in uniprot_kinases_short]
    for file in files:
        if any(file.startswith(uniprot_kinase) for uniprot_kinase in uniprot_kinases):
            if file.endswith('.annotations.txt.gz'):
                old_file_dir = str('/net/home.isilon/ds-russell/mechismoX/analysis/mutations/data/VLatest/'+ file[:4]+ '/'+ file)
                new_file_dir = str('/net/home.isilon/ds-russell/kinaseResistance/nat_vars/'+file)
                #print(old_file_dir)
                shutil.copyfile(old_file_dir, new_file_dir )
                print(file)


