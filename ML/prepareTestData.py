#!/usr/bin/env python
# coding: utf-8

import fetchData
import pickle
import pandas as pd
from cls import Kinase, Mutation

PTM_TYPES = ['ac', 'gl', 'm1', 'm2', 'm3', 'me', 'p', 'sm', 'ub']
AA = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
WS = 5

def predict(inputFile = 'sample_mutations.txt', BASE_DIR = './'):
    mydb = fetchData.connection()
    mydb.autocommit = True
    mycursor = mydb.cursor()

    trainMat = 'Acc\tGene\tMutation\tDataset\t'
    trainMat += 'hmmPos\thmmSS\thmmScoreWT\thmmScoreMUT\thmmScoreDiff\t'
    trainMat += 'Phosphomimic\t'
    trainMat += 'ChargesWT\tChargesMUT\tChargesDiff\t'
    startWS = int((WS-1)/2) * -1
    endWS = int((WS-1)/2)
    for position in range(startWS, endWS+1):
        trainMat += ('_'+str(position)+'\t').join(PTM_TYPES) + '_'+str(position)+'\t'
        trainMat += ('_'+str(position)+'_pfam\t').join(PTM_TYPES) + '_' + str(position) + '_pfam\t'
    trainMat += '_WT\t'.join(AA) + '_WT\t'
    trainMat += '_MUT\t'.join(AA) + '_MUT\t'
    trainMat += '\t'.join(['allHomologs','exclParalogs','specParalogs','orthologs','bpso','bpsh']) + '\t'
    for position in range(startWS, endWS+1):
        trainMat += ('_'+str(position)+'\t').join(['A', 'D', 'R']) + '_'+str(position)+'\t'
    trainMat += '\n'

    kinases = {}
    for line in open(inputFile, 'r'):
        if line[0] in ['#', '\n'] or len(line.split()) == 0: continue
        row = []
        # acc = line.split()[0]
        # gene = line.split()[1]
        print (line.split())
        name = line.split()[0].split('/')[0]
        acc, gene, uniprot_id = fetchData.getAccGene(mycursor, name)
        if acc not in kinases:
            kinases[acc] = Kinase(acc, gene)
        mutation = line.split('/')[1].rstrip()
        if mutation not in kinases[acc].mutations:
            kinases[acc].mutations[mutation] = Mutation(mutation, '-', acc, 'test')
        position = int(mutation[1:-1])
        mutAA = mutation[-1]
        wtAA = mutation[0]
        hmmPos, hmmScoreWT, hmmScoreMUT, hmmSS = fetchData.getHmmPkinaseScore(mycursor, acc, wtAA, position, mutAA)
        if hmmPos == '-': continue
        ptm_row = fetchData.getPTMscore(mycursor, acc, position, WS)
        aa_row = fetchData.getAAvector(wtAA, mutAA)
        homology_row = fetchData.getHomologyScores(mycursor, acc, wtAA, position, mutAA)
        if homology_row == None: continue
        is_phosphomimic = kinases[acc].mutations[mutation].checkPhosphomimic()
        charges_row = kinases[acc].mutations[mutation].findChangeInCharge()
        adr_row = fetchData.getADRvector(mycursor, acc, position, kinases, WS)
        ##save
        trainMat += acc + '\t' + gene + '\t' + mutation + '\t' + 'test' + '\t'
        trainMat += str(hmmPos) + '\t' + str(hmmSS) + '\t' + str(hmmScoreWT) + '\t' + str(hmmScoreMUT) + '\t' + str(hmmScoreMUT-hmmScoreWT) + '\t'
        trainMat += str(is_phosphomimic) + '\t'
        trainMat += '\t'.join([str(item) for item in charges_row]) + '\t'
        trainMat += '\t'.join([str(item) for item in ptm_row]) + '\t'
        trainMat += '\t'.join([str(item) for item in aa_row]) + '\t'
        trainMat += '\t'.join([str(item) for item in homology_row]) + '\t'
        trainMat += '\t'.join([str(item) for item in adr_row]) + '\n'

    with open(BASE_DIR+'/ML/'+'test_data.tsv', 'w') as f:
        f.write(trainMat)

    df = pd.read_csv(BASE_DIR+'/ML/'+'test_data.tsv', sep = '\t')
    # include columns
    columns_to_consider = []
    for line in open(BASE_DIR+'/ML/'+'columns_to_consider.txt', 'r'):
        columns_to_consider.append(line.strip())

    # df = df.loc[:, df.columns.isin(columns_to_consider)]
    test_data = df.to_numpy()
    features = df.loc[:, df.columns.isin(columns_to_consider)]
    features = features.to_numpy()
    features = features[:, 3:]

    filename = BASE_DIR+'/ML/'+'finalized_model.sav'
    clf = pickle.load(open(filename, 'rb'))
    scaler = pickle.load(open(BASE_DIR+'/ML/'+'finalized_scaler.pkl', 'rb'))
    features = scaler.transform(features)
    results = {}
    for row, predict in zip(test_data, clf.predict_proba(features)):
        print (row[0], row[1], row[2], row[4], round(predict[1], 3))
        name = row[0] + '/' + row[2]
        results[name] = {'prediction':round(predict[1], 3), 'hmmPos':row[4]}
    return results


if __name__ == '__main__':
    predict('sample_mutations.txt')