#/usr/bin/env python3
# coding: utf-8

'''
Calculate the log-odd scores
from a HMM file and an aln file
'''

import gzip, argparse
import numpy as np

def calculate_log_odds(input_hmm):
    '''
    A function that takes HMM file as input
    and returns 2 dics: log_odd scores and
    prevalent aa for an alignment position
    '''
    dic_background = {}
    dic_log_odds = {}
    dic_prevalent_aa = {}
    for line in gzip.open(input_hmm, 'rt'):
        if line.split()[0] == 'COMPO':
            background_scores = line.replace('\n', '').split()[1:]
            # print (background_scores)
            for background_score, aa in zip(background_scores, AA):
                dic_background[aa] = float(background_score)
        elif line.split()[0] == 'HMM':
            AA = line.replace('\n', '').split()[1:]
        elif line.split()[-1] == '-' and line.split()[-2] == '-':
            aln_position = int(line.split()[-5])
            dic_prevalent_aa[aln_position] = line.split()[-4]
            scores = line.replace('\n', '').split()[1:21]
            if aln_position not in dic_log_odds: dic_log_odds[aln_position] = {}
            for score, aa in zip(scores, AA):
                dic_log_odds[aln_position][aa] = round(float(score) - dic_background[aa], 5)
            # print (aln_position)
    return dic_log_odds, dic_prevalent_aa

class Protein:
    '''
    A class to store protein information
    '''
    def __init__(self, id) -> None:
        self.id = id
        self.seq = ''

def make_alignment(input_aln):
    '''
    A function that takes CLUSTAL-formatted alignment
    file as input and creates a dic that stores Protein
    class info. It also convert the alignent into a
    matrix format and returns the name of human id in
    the alignment
    '''
    dic_proteins = {}
    alignment = []
    human_name = ''
    for line in gzip.open(input_aln, 'rt'):
        if line.split() == []:
            continue
        if line[0] == '#':
            continue
        # print (line.split())
        if line.split()[0] in ['CLUSTAL', '//']:
            continue
        uniprot_acc = line.split()[0]
        if 'HUMAN_' in uniprot_acc and human_name == '': # HUMAN_* can occur more than once, esp in paralogs - take the first one
            human_name = uniprot_acc
        seq_aln = line.replace('\n', '').split()[1].upper()
        # print (seq_aln)
        if uniprot_acc not in dic_proteins: dic_proteins[uniprot_acc] = Protein(uniprot_acc)    
        dic_proteins[uniprot_acc].seq += seq_aln

    for uniprot_acc in dic_proteins:
        row = []
        for char in dic_proteins[uniprot_acc].seq:
            row.append(char)
        alignment.append(row)

    alignment = np.array(alignment)
    # print (alignment[:,9])
    return dic_proteins, alignment, human_name

def main(input_aln, input_hmm):
    '''
    The main function that calls the above
    2 functions and returns the output text
    Execute this when file runs as module
    and not when run as script
    '''
    dic_log_odds, dic_prevalent_aa = calculate_log_odds(input_hmm=input_hmm)
    dic_proteins, alignment, human_name = make_alignment(input_aln=input_aln)

    count = 0; text = ''
    for aln_position, wt_aa in enumerate(dic_proteins[human_name].seq):
        if wt_aa in ['-', '.']:
            continue
        if aln_position+1 not in dic_log_odds:
            print ('aln_position', aln_position, 'does not exist in profile HMM')
            continue
        count += 1
        # print (aln_position+1, wt_aa, count)
        # print (dic_log_odds.keys())
        wt_score = dic_log_odds[aln_position+1][wt_aa]
        for mut_aa in 'ACDEFGHIKLMNPQRSTVWY':
            mut_score = dic_log_odds[aln_position+1][mut_aa]
            text += human_name.split('_')[1] + '/' + wt_aa + str(count) + mut_aa + ' '
            text += dic_prevalent_aa[aln_position+1] + ' '
            text += str(wt_score) + ' ' + str(mut_score) + ' ' + str(round(mut_score-wt_score, 5)) + ' '
            text += ''.join(alignment[:,aln_position])  + '\n'

    print(text)
    return text

if __name__ == '__main__':
    '''
    Execute this when file runs as script
    and not when run as module
    '''
    parser = argparse.ArgumentParser(description='Generate log-odds score from *.aln.gz and corresponding *.hmm.gz files',
                                        epilog='gurdeep.singh[at]bioquant[.]uni[-]heidelberg[.]de',
                                        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('aln', help='Alignment file in CLUSTAL W/STOCKHOLM formats')
    parser.add_argument('hmm', help='profile HMM file')
    args = parser.parse_args()
    input_aln = args.aln
    input_hmm = args.hmm
    main(input_aln, input_hmm)
