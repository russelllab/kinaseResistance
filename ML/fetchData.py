#!/usr/bin/env python
# coding: utf-8

import os, gzip

'''
List of functions that fetch data from
different files
'''

def fetchFasta(kinases, Kinase):
    for line in open('../data/humanKinases.fasta', 'r'):
        #print (line)
        if line[0] == '>':
            acc = line.split('|')[1].replace('\n', '')
            gene = line.split('GN=')[1].split()[0]
            kinases[acc] = Kinase(acc, gene)
            # flag = 0
            # if acc not in exceptions:
            #     kinases[acc] = kinase(acc, gene)
            #     flag = 1
        else:
            # if flag == 1:
            kinases[acc].fasta += line.replace('\n', '')

def fetchGroup(kinases, Kinase):
    for line in open('../data/kinases.tsv', 'r'):
        acc = line.split('\t')[7].split('>')[1].split('<')[0]
        if acc in kinases:
            kinases[acc].group = line.split('\t')[4]

def fetchPkinaseHMM():
    hmm = {} # hmmPosition > AA > bit-score
    for line in open('../pfam/Pkinase.hmm'):
        if len(line.split()) > 2:
            if line.split()[-2] == '-' and line.split()[-3] == '-':
                #print (line.split())
                position = int(line.split()[0])
                hmm[position] = {}
                for value, aa in zip(line.split()[1:-5], AA):
                    hmm[position][aa] = float(value)
            elif line.split()[0] == 'HMM':
                AA = line.replace('\n', '').split()[1:]
    return hmm

def fetchHmmsearch(kinases, Kinase):
    '''
    Function to do an hmmsearch of all kinases against Pkinase.hmm
    and store the mappings. Note that some kinases may have more than
    one Pkinase domain.
    '''
    os.system('hmmsearch -o out.txt ../pfam/Pkinase.hmm ../data/humanKinases.fasta')
    flag = 0
    for line in open('out.txt', 'r'):
        if line[:2] == '>>':
            acc = line.split('|')[1]
            flag = 1
            #print (acc)
        if flag == 1 and line.split()!= [] and acc in kinases:
            if '== domain' in line:
                domainNum = int(line.split('domain')[1].split()[0])
                kinases[acc].domains[domainNum] = {}
            elif line.split()[0] == 'Pkinase':
                hmmStart = int(line.split()[1])
                hmmSeq = line.split()[2]
                hmmEnd = int(line.split()[3])
            elif acc in line.split()[0]:
                kinaseStart = int(line.split()[1])
                kinaseSeq = line.split()[2]
                kinaseEnd = int(line.split()[3])
                for hmmChar, kinaseChar in zip(hmmSeq, kinaseSeq):
                    if hmmChar not in ['.', '-'] and kinaseChar not in ['.', '-']:
                        #kinases[acc].domains[domainNum][kinaseStart] = hmmStart
                        kinases[acc].domains[domainNum][hmmStart] = kinaseStart
                        hmmStart += 1
                        kinaseStart += 1
                    elif hmmChar in ['.', '-']:
                        kinaseStart += 1
                    elif kinaseChar in ['.', '-']:
                        hmmStart += 1
        #print (kinases[acc].domains)
        #sys.exit()
    print (kinases['Q96NX5'].domains[1][174])

def createDicForDSSP(dic, position, mutation, value):
    if position not in dic:
        dic[position] = {}
    dic[position][mutation] = float(value)

def dsspScores(kinases, Kinase):
    remove = []
    dir = '/net/home.isilon/ds-russell/mechismoX/analysis/features/data/VLatest/'
    for num, acc in enumerate(kinases):
        for line in gzip.open(dir + acc[:4] + '/AF-' + acc + '-F1-model_v1.dssp-scores.gz', 'rt'):
            if line.split()[0] == 'length':
                continue
            #print (acc, line.split())
            position = int(line.split()[0])
            mutation = line.split()[2] + line.split()[0] + line.split()[10]
            ## Dihedral angles
            torsional = line.split()[18]
            #print (kinases[acc].dihedral)
            createDicForDSSP(kinases[acc].dihedral, position, mutation, torsional)
            ## Secondary structures
            secondary = line.split()[22]
            createDicForDSSP(kinases[acc].sec, position, mutation, secondary)
            ## Accessibility
            accessibility = line.split()[26]
            createDicForDSSP(kinases[acc].access, position, mutation, accessibility)
            ## Buried
            buried = line.split()[30]
            createDicForDSSP(kinases[acc].burr, position, mutation, buried)
            #break
        #break
    print (list(set(remove)))

def iupredScores(kinases, Kinase):
    remove = []
    dir = '/net/home.isilon/ds-russell/mechismoX/analysis/features/data/VLatest/'
    for num, acc in enumerate(kinases):
        if ((num+1)%50 == 0):
            print (num+1)
        for line in gzip.open(dir + acc[:4] + '/AF-' + acc + '-F1-model_v1.iupred.gz', 'rt'):
            #print (acc, line.split())
            position = int(line.split()[0])
            mutation = line.split()[2] + line.split()[0] + line.split()[3]
            ## IUPred
            iupred = float(line.split()[9])
            createDicForDSSP(kinases[acc].iupred, position, mutation, iupred)


def mechismoScores(kinases, Kinase):
    remove = []
    dir = '/net/home.isilon/ds-russell/mechismoX/analysis/features/data/VLatest/'
    for num, acc in enumerate(kinases):
        if ((num+1)%50 == 0):
            print (num+1)
        for line in gzip.open(dir + acc[:4] + '/AF-' + acc + '-F1-model_v1.mech_intra.gz', 'rt'):
            if line.split()[0] == 'MECH':
                #print (acc, line.split())
                #sys.exit()
                position = int(line.split()[1])
                mutation = line.split()[2] + line.split()[1] + line.split()[3]
                ## Mechismo score
                mechismo = float(line.split()[6])
                createDicForDSSP(kinases[acc].mechismo, position, mutation, mechismo)

def homologyScores(kinases, Kinase):
    remove = []
    path = '/net/home.isilon/ds-russell/mechismoX/analysis/alignments/data/HUMAN/orthologs_only/'
    for num, acc in enumerate(kinases):
        if ((num+1)%50 == 0):
            print (num+1)
        for dic, fileEnd in zip([
                            kinases[acc].allHomologs,
                            kinases[acc].orthologs,
                            kinases[acc].exclParalogs,
                            kinases[acc].specParalogs,
                            kinases[acc].bpso,
                            kinases[acc].bpsh],
                            [
                            '_all_homs.scores.txt.gz',
                             '_orth.scores.txt.gz',
                             '_excl_para.scores.txt.gz',
                             '_spec_para.scores.txt.gz',
                             '_bpso.scores.txt.gz',
                             '_bpsh.scores.txt.gz'
                             ]):
            if os.path.isfile(path + acc[:4] + '/' + acc + fileEnd) is False:
                print (path + acc[:4] + '/' + acc + fileEnd, 'does not exist')
                continue
            for line in gzip.open(path + acc[:4] + '/' + acc + fileEnd, 'rt'):
                #print (acc, line.split())
                #sys.exit()
                value = line.split()[0].split('/')[1]
                position = value[1:-1]
                residue = value[-1]
                #print (mutation, position)
                ## Mechismo score
                score = float(line.split()[4])
                createDicForDSSP(dic, position, residue, score)
