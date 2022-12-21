#!/usr/bin/env python3
# coding: utf-8
'''
Script to perform a test input mutations
'''
import os, sys, argparse, gzip, math
# from turtle import position
import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact
# import plotly.express as px

class Mutation():
    def __init__(self, mutation):
        self.mutation = mutation
        self.scores = {}

class Kinase():
    def __init__(self, name):
        self.name = name
        self.acc = ''
        self.gene = ''
        self.seq = ''
        self.ptmsites = {}
        self.act_deact_res = {}
        self.mutations = {}
        self.hmm2seq = {}
        self.fasta2pfam = {}
        self.fetch_fasta_sequence()
        self.fetch_hmmsearch()
        self.fetch_acc_gene()
    
    def fetch_acc_gene(self):
        pass

    def fetch_fasta_sequence(self):
        '''
        Fetch FASTA sequence of the kinase
        '''
        for line in gzip.open('../KA/UniProtFasta2/'+self.name+'.fasta.gz', 'rt'):
            if line[0] != '>':
                self.seq += line.replace('\n', '')
            else:
                self.gene = line.split('GN=')[1].split()[0]
    
    def fetch_hmmsearch(self):
        '''
        Fetch HMMsearch of the kinase
        '''
        flag = 0
        for line in open('humanKinasesHmmsearch.txt', 'r'):
            if line[:2] == '>>':
                if flag == 0: flag = 1
                if flag == 2: break
                if self.name in line: flag += 1
                continue
            if flag != 2:
                continue
            if len(line.split()) == 0:
                continue
            
            # print (self.name.split(), line.split()[0])
            if 'Pkinase' in line:
                pfam_start, pfam_seq, pfam_end = int(line.split()[1]), line.split()[2], int(line.split()[3])
            elif self.name in line:
                fasta_start, fasta_seq, fasta_end = int(line.split()[1]), line.split()[2], int(line.split()[3])
                dic_pfam2fasta = {}
                dic_fasta2pfam = {}
                for pfam_aa, fasta_aa in zip(pfam_seq, fasta_seq):
                    if pfam_aa not in ['.', '-'] and fasta_aa not in ['.', '-']:
                        self.fasta2pfam[fasta_start] = pfam_start
                        pfam_start += 1
                        fasta_start += 1
                    elif pfam_aa not in ['.', '-']:
                        pfam_start += 1
                    elif fasta_aa not in ['.', '-']:
                        fasta_start += 1
        # print (self.fasta2pfam[315])
        # sys.exit()

def get_ptm_sites(dic_kinases):
    '''
    Fetch PTM sites
    '''
    dic_ptm_sites = {}
    for line in open('Kinase_psites4.tsv', 'r'):
        if line[0] == '#':
            continue
        if line.split('\t')[2] != 'Pkinase':
            continue
        position = int(line.split('\t')[4])
        acc = line.split('\t')[0]
        ptm_site = int((line.split('\t')[3].split('-')[0])[1:])
        ptm_type = line.split('\t')[3].split('-')[1]
        if position not in dic_ptm_sites: dic_ptm_sites[position] = {}
        if ptm_type not in dic_ptm_sites[position]: dic_ptm_sites[position][ptm_type] = {}
        if acc not in dic_ptm_sites[position][ptm_type]: dic_ptm_sites[position][ptm_type][acc] = []
        dic_ptm_sites[position][ptm_type][acc].append(ptm_site)
        if acc in dic_kinases:
            # print (acc, line, dic_kinases[acc].ptmsites)
            if ptm_type not in dic_kinases[acc].ptmsites:
                dic_kinases[acc].ptmsites[ptm_type] = []
            dic_kinases[acc].ptmsites[ptm_type].append(ptm_site)
    
    # print (dic_ptm_sites[156]['p'])
    return dic_ptm_sites

def get_act_deact_res_sites(dic_kinases):
    '''
    Fetch act/deact/res sites
    '''
    dic_act_deact_res_sites = {}
    for line in open('../KA/act_deact_mut_for_scores_fin.tsv', 'r'):
        if line.split('\t')[0] == 'uniprot_id':
            continue
        if 'de' in line.split('\t')[2] or 'dup' in line.split('\t')[2]:
            continue
        acc = line.split('\t')[0]
        position = int(line.split('\t')[2])
        wt_aa = line.split('\t')[1]
        mut_aa = line.split('\t')[3]
        mut_type = line.split('\t')[-1].replace('\n', '')
        if acc not in dic_act_deact_res_sites: dic_act_deact_res_sites[acc] = {'R':[], 'A':[], 'D':[]}
        dic_act_deact_res_sites[acc][mut_type].append(position)
        if acc in dic_kinases:
            # print (acc, line, dic_kinases[acc].ptmsites)
            if mut_type not in dic_kinases[acc].act_deact_res: dic_kinases[acc].act_deact_res[mut_type] = {}
            if position not in dic_kinases[acc].act_deact_res[mut_type]: dic_kinases[acc].act_deact_res[mut_type][position] = []
            dic_kinases[acc].act_deact_res[mut_type][position].append(wt_aa)
    
    for line in open('../KA/resistant_mutations_Nov22.tsv.gz', 'r'):
        if line[0] == '#':
            continue
        if 'del' in line.split('\t')[2] or 'du' in line.split('\t')[2] or '_' in line.split('\t')[2]:
            continue
        acc = line.split('\t')[1]
        position = int((line.split('\t')[2])[1:-1])
        wt_aa = (line.split('\t')[2])[0]
        mut_aa = (line.split('\t')[2])[0]
        mut_type = 'R'
        if acc not in dic_act_deact_res_sites: dic_act_deact_res_sites[acc] = {'R':[], 'A':[], 'D':[]}
        dic_act_deact_res_sites[acc][mut_type].append(position)
        if acc in dic_kinases:
            # print (acc, line, dic_kinases[acc].ptmsites)
            if mut_type not in dic_kinases[acc].act_deact_res: dic_kinases[acc].act_deact_res[mut_type] = {}
            if position not in dic_kinases[acc].act_deact_res[mut_type]: dic_kinases[acc].act_deact_res[mut_type][position] = []
            dic_kinases[acc].act_deact_res[mut_type][position].append(wt_aa)
    
    # print (dic_ptm_sites[156]['p'])
    return dic_act_deact_res_sites

def get_log_odds(dic_kinases, alignment_types):
    '''
    '''
    # alignment_types = ['all_homs','bpsh','bpso','excl_para','orth','spec_para']
    path_to_alignments = '/net/home.isilon/ds-russell/mechismoX/analysis/alignments/data/HUMAN/orthologs_only/'
    for name in dic_kinases:
        for aln_type in alignment_types:
            if os.path.isfile(path_to_alignments+name[:4]+'/'+name+'_'+aln_type+'.scores.txt.gz') is False:
                continue
            for line in gzip.open(path_to_alignments+name[:4]+'/'+name+'_'+aln_type+'.scores.txt.gz', 'rt'):
                if line[0] == '\n':
                    continue
                mutation = line.split()[0].split('/')[1]
                if mutation in dic_kinases[name].mutations:
                    object_mutation = dic_kinases[name].mutations[mutation]
                    object_mutation.scores[aln_type] = float(line.split()[4])
            # print (mutation)

def do_BLAST(name, pfam_pos, ptm_type, dic_ptm_sites):
    # print (pfam_pos, ptm_type, dic_ptm_sites)
    seq = ''
    for acc in dic_ptm_sites[pfam_pos][ptm_type]:
        # print (acc)
        for line in gzip.open('../KA/UniProtFasta2/'+acc+'.fasta.gz', 'rt'):
            if line[0] == '>':
                seq += '>' + line.split('GN=')[1].split()[0] + '\n'
            else:
                seq += line
    # print (seq)
    open('blastdb/inputDB.fasta', 'w').write(seq)
    os.system('makeblastdb -in blastdb/inputDB.fasta -dbtype prot -out blastdb/inputDB  -logfile blastdb/inputDB.log')
    os.system('gunzip ../KA/UniProtFasta2/'+name+'.fasta.gz')
    os.system('blastp -query ../KA/UniProtFasta2/'+name+'.fasta -db blastdb/inputDB -outfmt=6 -out outputDB.txt')
    os.system('gzip ../KA/UniProtFasta2/'+name+'.fasta')
    for line in open('outputDB.txt', 'r'):
        best_hit_gene = line.split('\t')[1]
        best_hit_identity = str(math.floor(float(line.split('\t')[2])))
        break
    # return best_hit.split('|')[2]
    return best_hit_gene, best_hit_identity

def main(dic_kinases):
    '''
    '''
    dic_ptm_sites = get_ptm_sites(dic_kinases)
    dic_act_deact_res_sites = get_act_deact_res_sites(dic_kinases)
    alignment_types = ['all_homs','bpsh','bpso','excl_para','orth','spec_para']
    get_log_odds(dic_kinases, alignment_types)
    # print (dic_kinases['Q02750'].ptmsites)
    # print (dic_ptm_sites)
    # prepare_fasta_sequence(dic_kinases)
    ptm_types = ['ac','gl','m1','m2','m3','me','p','sm','ub']
    mut_types = ['A', 'D', 'R']
    heading = '#Gene\tName\tMutation\tPfam_Pos\t'
    heading += '\t'.join(ptm_types) + '\t'
    heading += '\t'.join(ptm_types) + '\t'
    heading += 'Pfam_PTM' + '\t'
    heading += '\t'.join(mut_types) + '\t'
    heading += '\t'.join(alignment_types) + '\t'
    heading += 'MechismoX URL'
    print (heading)
    for name in dic_kinases:
        for mutation in dic_kinases[name].mutations:
            # print (name, mutation, dic_kinases[name].mutations[mutation].ortho_scores)
            # sys.exit()
            wt_aa = mutation[0]
            mut_aa = mutation[-1]
            fasta_pos = int(mutation[1:-1])
            row = []
            ## Ignore if the given position CANNOT be mapped to Pkinase
            if fasta_pos not in dic_kinases[name].fasta2pfam:
                continue
            ## Check if +1/-1 residues are PTM-sites
            # for ptm_type in dic_kinases[name].ptmsites:
            for ptm_type in ptm_types:
                if ptm_type not in dic_kinases[name].ptmsites:
                    row.append('-')
                    continue
                if fasta_pos+1 in dic_kinases[name].ptmsites[ptm_type]:
                    # print (name, dic_kinases[name].gene, mutation, '-', fasta_pos+1, 'is a', ptm_type+'-site')
                    row.append('Y')
                    continue
                if fasta_pos-1 in dic_kinases[name].ptmsites[ptm_type]:
                    # print (name, dic_kinases[name].gene, mutation, '-', fasta_pos-1, 'is a', ptm_type+'-site')
                    row.append('Y')
                    continue
                row.append('-')
            ## Check if it is a PTM site in the given protein
            pfam_pos = dic_kinases[name].fasta2pfam[fasta_pos]
            # for ptm_type in dic_ptm_sites[pfam_pos]:
            for ptm_type in ptm_types:    
                if pfam_pos in dic_ptm_sites:
                    if ptm_type not in dic_ptm_sites[pfam_pos]:
                        row.append('-')
                        continue
                    if name not in dic_ptm_sites[pfam_pos][ptm_type]:
                        row.append('-')
                        continue
                    # print (pfam_pos, acc, name, dic_ptm_sites[pfam_pos][ptm_type][acc])
                    if fasta_pos not in dic_ptm_sites[pfam_pos][ptm_type][name]:
                        row.append('-')
                        continue
                    row.append('Y')
                    # print (name, dic_kinases[name].gene, mutation, pfam_pos, fasta_pos, 'is a', ptm_type+'-site')
                else:
                    row.append('-')
            ## Check if it is a PTM site at the given PFAM position
            ptms = []
            if pfam_pos in dic_ptm_sites:
                for ptm_type in ptm_types:
                    if ptm_type in dic_ptm_sites[pfam_pos]:
                        best_hit_gene, best_hit_identity = do_BLAST(name, pfam_pos, ptm_type, dic_ptm_sites)
                        ptms.append(best_hit_gene+'|'+ptm_type+'|'+best_hit_identity)
                        # print (name)
            if (len(ptms) > 0):
                row.append(';'.join(ptms))
            else:
                row.append('-')
            ## Act/Deact/Res sites
            # for mut_type in dic_kinases[name].act_deact_res:
            for mut_type in mut_types:
                if mut_type not in dic_kinases[name].act_deact_res:
                    row.append('-')
                    continue
                if fasta_pos not in dic_kinases[name].act_deact_res[mut_type]:
                    row.append('-')
                    continue
                row.append(mut_type)
                # print (name, dic_kinases[name].gene, mutation, pfam_pos, fasta_pos, 'is a/an', mut_type+'-site')
            ## Check if there are scores
            for aln_type in alignment_types:
                if aln_type in dic_kinases[name].mutations[mutation].scores:
                    row.append(str(dic_kinases[name].mutations[mutation].scores[aln_type]))
                else:
                    row.append('-')
            row.append('http://mechismox.russelllab.org/result?protein='+name+'&mutation='+mutation)
            
            print (dic_kinases[name].gene + '\t' + name + '\t' + mutation + '\t' + str(pfam_pos) + '\t' + '\t'.join(row))
            # if len(row) == 22: sys.exit()

if __name__ == '__main__':
    '''
    Execute this when file runs as script
    and not when run as module
    '''
    parser = argparse.ArgumentParser(description='Test Kinase mutations',
                                        epilog='gurdeep.singh[at]bioquant[.]uni[-]heidelberg[.]de',
                                        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('input', help='File with input in Mechismo format')
    # parser.add_argument('hmm', help='profile HMM file')
    args = parser.parse_args()
    # input_muts = args.input
    dic_kinases = {}
    for line in open(args.input, 'r'):
        if line[0] == '#':
            continue
        name = line.split('/')[0]
        mutation = line.split('/')[1].replace('\n', '')
        if name not in dic_kinases: dic_kinases[name] = Kinase(name)
        dic_kinases[name].mutations[mutation] = Mutation(mutation)
    # input_hmm = args.hmm
    main(dic_kinases)
