import os
import pandas as pd
from Bio import SeqIO
import gzip
#import shutil
#select from robScores2 the kinase uniprot IDs and save them as a list for parcing in two formats once as
#first 4 letters and another as full
uniprot_kinases_raw = os.listdir('/net/home.isilon/ds-russell/kinaseResistance/KA/robScores2/')
#print(uniprot_kinases_raw)
uniprot_kinases = []
uniprot_kinases_short = []
for uniprot_kinase in uniprot_kinases_raw:
    if uniprot_kinase.endswith('.tsv.gz'):
        uniprot_kinase_id = uniprot_kinase.split('.')[0]
        uniprot_kinases.append(uniprot_kinase_id)
        uniprot_kinase_id_short = uniprot_kinase[:4]
        uniprot_kinases_short.append(uniprot_kinase_id_short)
  
#create a dictionary with uniprot_type:associated_sequences #for ortho
orth_lst = {}
orth_cutoff = {}
for root, dirs, files in os.walk('/net/home.isilon/ds-russell/kinaseResistance/kin_seq/'):
    for file in files:
        if file.endswith('_orth.aln.gz') is False:
            continue
        #print(file)
        df = pd.read_csv(str('/net/home.isilon/ds-russell/kinaseResistance/kin_seq/'+file),
                                   sep = '\s+|\t+|\s+\t+|\t+\s+', engine = 'python')
        #print(df)
        df['uniprot_id'] = df.iloc[:,0]
        df = df['uniprot_id'].str.strip()
        #add first row
        if df[0] != str('HUMAN_' + file.split('_')[0]):
            first_value = pd.Series([str('HUMAN_' + file.split('_')[0])])
            #print(first_value)
            df = pd.concat([first_value, df], axis = 0, ignore_index = True)


        df = df.drop_duplicates()
        uniprot_type = str(file.split('.')[0])
        
        #org_unirot associates ids 
        orth_lst[uniprot_type] = list(df)
        
        #cutoff 
        tfp_nrows = round(len(df) * 0.25) #25%
        fp_nrows = round(len(df) * 0.50) #50%
        sfp_nrows = round(len(df) * 0.75) #75%
        hp_rows = len(df) #100%
        cutoff = [tfp_nrows,fp_nrows,sfp_nrows, hp_rows]
        orth_cutoff[uniprot_type] = cutoff
        #print(len(df[:20]))
        #print(file)
        
            
#loop over fa.gz files  and create a dictionary with organism_uniprot:sequence
dict_glob_seq = {} #for all files

for root, dirs, files in os.walk('/net/home.isilon/ds-russell/kinaseResistance/kin_seq/'):
    #all scope of files
    for file in files:
        #within a file
        if file.endswith('fa.gz') is False:
            continue
        #print(file)
        with gzip.open("/net/home.isilon/ds-russell/kinaseResistance/kin_seq/"+file, "rt") as handle:
            
            lines = handle.readlines() #saves each line as an item in a list
            #print(len(lines))
            dict_fasta = {} #for each file
            
            for line in lines:
                if line.startswith('>'):
                    name = line.split('>')[1].replace('\n', '')
                    #print (name)
                    dict_fasta[name] = ''
                    
                else:
                    #print (line)
                    dict_fasta[name] += line.replace('-', '').replace('\n', '')
                    
                    
            dict_glob_seq.update(dict_fasta)
            
            
#25 % ortho
for uniprotid_type, assoc_uniprotids in orth_lst.items():
    print(uniprotid_type)
    fasta = ''
    file_w_path = '/net/home.isilon/ds-russell/kinaseResistance/AK_output/'+ str(uniprotid_type.split('.')[0]) + "_25_seq.fa.gz"
    
    for assoc_uniprotid in assoc_uniprotids[:orth_cutoff[uniprotid_type][0]]: 
        fasta += '>' + assoc_uniprotid + '\n' + dict_glob_seq[assoc_uniprotid]+ '\n'
    
    gzip.open(file_w_path, 'wt').write(fasta)
    