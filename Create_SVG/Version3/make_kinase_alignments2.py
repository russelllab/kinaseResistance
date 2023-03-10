'''
A script to make kinase alignments
It will take a kinase as input and 
do a BLAST against other kinases.
It will then pick the top n candidates
based on number of ADR mutations around
the position of interest. And finally
return a multiple sequence alignment.
'''

import os, sys, gzip, pickle
import pandas as pd
PATH_TO_FASTA = '../../KA/UniProtFasta2/'
kinases = {}

class Kinase:
    '''
    Class to store kinase information
    '''
    def __init__(self, acc, gene):
        self.acc = acc
        self.gene = gene
        self.sequence = ''
        self.mut_types = {'A':[], 'R': [], 'D':[]}
        self.seq2pfam = {}
        self.pfam2seq = {}

def loadKinaseInfo():
    '''
    Function load kinase sequence and information in the kinases dic (class)
    '''
    ## FASTA
    for files in os.listdir(PATH_TO_FASTA):
        # ignore files that do not end with .fasta.gz
        if files.endswith('.fasta.gz') is False:
            continue
        # Ignore isoforms
        if '-' in files:
            continue
        for line in gzip.open(PATH_TO_FASTA+files, 'rt'):
            if line[0] == '>':
                acc = line.split('|')[1]
                gene = line.split('GN=')[1].split()[0]
                if acc not in kinases:
                    kinases[acc] = Kinase(acc, gene)
            else:
                kinases[acc].sequence += line.replace('\n', '')
    
    ## Act/deact Information
    for line in open('../../AK_mut_w_sc_feb2023/act_deact_v2.tsv', 'r'):
        if line.split('\t')[0] == 'uniprot_name': continue
        acc = line.split('\t')[1]
        position = str(line.split('\t')[3])
        mut_type = line.split('\t')[5]
        if position not in kinases[acc].mut_types[mut_type]:
            kinases[acc].mut_types[mut_type].append(position)
    
    ## Resistant Information
    for line in gzip.open('../../KA/resistant_mutations_Mar_2023.tsv.gz', 'rt'):
        if line[0] == '#': continue
        acc = line.split('\t')[2]
        position = str(line.split('\t')[5].replace('\n', ''))
        mut_type = 'R'
        if position not in kinases[acc].mut_types[mut_type]:
            kinases[acc].mut_types[mut_type].append(position)
    
    ## Pfam2Seq mappings
    for line in gzip.open('../../data/allKinasesPkinaseMappings.tsv.gz', 'rt'):
        if line[0] == '#': continue
        acc = line.split('\t')[0]
        if '-' in acc: continue
        if acc not in kinases: continue
        pfam_pos = line.split('\t')[4].replace('\n', '')
        seq_pos = line.split('\t')[2].replace('\n', '')
        kinases[acc].seq2pfam[seq_pos] = pfam_pos
        kinases[acc].pfam2seq[pfam_pos] = seq_pos
    
    # print (kinases['P46734'].mut_types)

def runBlast(input_kinase, n_candidates = 10, input_position = 174, window_size = 10):
    '''
    Function to run BLASTp against all kinases
    '''
    input_kinase_acc = None
    input_kinase_gene = None
    all_kinases = ''
    for acc in kinases:
        if acc == input_kinase or kinases[acc].gene == input_kinase:
            input_kinase_acc = acc
            input_kinase_gene = acc
        else:
            all_kinases += '>'+acc+'|'+kinases[acc].gene+'\n'
            all_kinases += kinases[acc].sequence + '\n'
    if input_kinase_acc == None or input_kinase_gene == None:
        print ('input kinase not found')
        sys.exit()

    # print (all_kinases)
    # Save all kinases sequence in FASTA format
    gzip.open('../BLASTdb/all_kinases.fasta.gz', 'wt').write(all_kinases)
    # Run makeblastdb
    os.system('gunzip ../BLASTdb/all_kinases.fasta.gz')
    os.system('makeblastdb -dbtype prot -in ../BLASTdb/all_kinases.fasta -out ../BLASTdb/all_kinases_db')
    os.system('gzip ../BLASTdb/all_kinases.fasta')
    # Run BLASTp
    os.system('gunzip '+PATH_TO_FASTA+input_kinase_acc+'.fasta.gz')
    os.system('blastp -out ../BLASTdb/output_blast.txt -outfmt 7 -query '+PATH_TO_FASTA+input_kinase_acc+'.fasta' + " -db ../BLASTdb/all_kinases_db")
    os.system('gzip '+PATH_TO_FASTA+input_kinase_acc+'.fasta')
    # Fetch top N candidates
    candidates = []
    for line in open('../BLASTdb/output_blast.txt', 'r'):
        if line[0]=='#': continue
        candidates.append(line.split('\t')[1])
        # if len(candidates) == n_candidates: break
    
    ## Eliminate cases where sum(ADR) mutations < 5
    ## at position 84 (say)
    # print (candidates)
    data = []
    start_window = int(window_size/2)
    end_window = int(window_size/2) + 1
    for num, candidate in enumerate(candidates, start=1):
        sum_ADR = 0
        row = []
        for mut_type in ['A', 'D', 'R']:
            for position in range(input_position-start_window, input_position+end_window):
                acc = candidate.split('|')[0]
                if str(position) in kinases[acc].mut_types[mut_type]:
                    sum_ADR += 1
        if sum_ADR >= 1:
            print (candidate, sum_ADR)
        row.append(candidate)
        row.append(sum_ADR)
        row.append(num)
        data.append(row)
    df = pd.DataFrame(data, columns = ['candidate', 'sum_ADR', 'identity'])
    df = df.sort_values(by=['sum_ADR', 'identity'], ascending=[False, False])
    print (df.head(n_candidates))
    candidates = df['candidate'].head(n_candidates).to_numpy()
    candidates = list(candidates)
    # print (candidates)
    # sys.exit()
    runHmmAlign(input_kinase, input_kinase_acc, input_kinase_gene, candidates)
    return input_kinase, input_kinase_acc, input_kinase_gene, candidates

def createDic(input_kinase, input_kinase_acc, input_kinase_gene, candidates):
    '''
    Function to reformat the hmmalign output to CLUSTAL
    and create a dictionary with mutation data and pfam
    regions data
    '''
    # reformat hmmalign output
    new_format = 'CLUSTAL\n\n'
    for line in open('../BLASTdb/output_hmmalign.aln', 'r'):
        if line[0] == '#': continue
        # if len(line.split()) == 0: continue
        if line[:2] == '//': continue
        new_format += line
    open(input_kinase + '.aln', 'w').write(new_format)
    
    # create dictionary for mutation data
    dic = {}
    for candidate in [input_kinase_acc]+candidates:
        acc = candidate.split('|')[0]
        if acc not in dic: dic[acc] = {}
        for mut_type in kinases[acc].mut_types:
            dic[acc][mut_type] = kinases[acc].mut_types[mut_type]
    with open(input_kinase+'_mutations.pkl', 'wb') as f:
        pickle.dump(dic, f)
    # with open(input_kinase+'_mutations.pkl', 'rb') as f:
    #     loaded_dic = pickle.load(f)
    # print (loaded_dic)

    # create dictionary for pfam regions data
    dic = {}
    for line in open('../../data/kinase_pfam_regions.tsv', 'r'):
        if line[0] == '#': continue
        region = line.split()[0]
        # print (line.split())
        start_pfam = line.split()[1]
        end_pfam = line.split()[2].replace('\n', '')
        # print (start_pfam, end_pfam, kinases[input_kinase_acc].pfam2seq)
        start_seq = kinases[input_kinase_acc].pfam2seq[start_pfam]
        end_seq = kinases[input_kinase_acc].pfam2seq[end_pfam]
        dic[region] = [start_seq, end_seq]
    with open(input_kinase+'_regions.pkl', 'wb') as f:
        pickle.dump(dic, f)
    # with open(input_kinase+'_regions.pkl', 'rb') as f:
    #     loaded_dic = pickle.load(f)
    # print (loaded_dic)

def runHmmAlign(input_kinase, input_kinase_acc, input_kinase_gene, candidates):
    # print(candidates)
    all_candidates = ''
    for candidate in [input_kinase_acc]+candidates:
        print (candidate)
        acc = candidate.split('|')[0]
        all_candidates += '>'+acc+'|'+kinases[acc].gene+'\n'
        all_candidates += kinases[acc].sequence + '\n'
    # Save all kinases sequence in FASTA format
    gzip.open('../BLASTdb/all_candidates.fasta.gz', 'wt').write(all_candidates)
    # Run BLASTp
    os.system('gunzip ../BLASTdb/all_candidates.fasta.gz')
    os.system('hmmalign -o ../BLASTdb/output_hmmalign.aln ../../pfam/Pkinase.hmm ../BLASTdb/all_candidates.fasta')
    os.system('gzip ../BLASTdb/all_candidates.fasta')
    createDic(input_kinase, input_kinase_acc, input_kinase_gene, candidates)
    return None

if __name__ == '__main__':
    input_kinase = 'PLK3'
    loadKinaseInfo()
    runBlast(input_kinase, n_candidates = 35, input_position=203, window_size = 31)