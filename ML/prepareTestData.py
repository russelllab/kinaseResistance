#!/usr/bin/env python3.10
# coding: utf-8

import fetchData
from tqdm import tqdm
import pickle
import pandas as pd
from cls import Kinase, Mutation
import argparse

PTM_TYPES = ['ac', 'gl', 'm1', 'm2', 'm3', 'me', 'p', 'sm', 'ub']
AA = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
WS = 5

def predict(inputFile, outputFile = None, BASE_DIR = '../'):
    '''
    Function to predict the effect of
    mutations on kinase activity
    '''

    # Connect to DB
    mydb = fetchData.connection()
    mydb.autocommit = True
    mycursor = mydb.cursor()

    # Make header
    data = []
    # trainMat = 'Acc\tGene\tMutation\tDataset\t'
    row = ['Acc','Gene','Mutation','Dataset']
    # trainMat += 'hmmPos\thmmSS\thmmScoreWT\thmmScoreMUT\thmmScoreDiff\t'
    row += ['hmmPos','hmmSS','hmmScoreWT','hmmScoreMUT','hmmScoreDiff']
    # trainMat += 'Phosphomimic\t'
    row += ['Phosphomimic']
    # trainMat += 'ChargesWT\tChargesMUT\tChargesDiff\t'
    row += ['ChargesWT','ChargesMUT','ChargesDiff']
    startWS = int((WS-1)/2) * -1
    endWS = int((WS-1)/2)
    for position in range(startWS, endWS+1):
        # trainMat += ('_'+str(position)+'\t').join(PTM_TYPES) + '_'+str(position)+'\t'
        row += [ptm_type+'_'+str(position) for ptm_type in PTM_TYPES]
        # trainMat += ('_'+str(position)+'_pfam\t').join(PTM_TYPES) + '_' + str(position) + '_pfam\t'
        row += [ptm_type+'_'+str(position)+'_pfam' for ptm_type in PTM_TYPES]
    # trainMat += '_WT\t'.join(AA) + '_WT\t'
    row += [aa+'_WT' for aa in AA]
    # trainMat += '_MUT\t'.join(AA) + '_MUT\t'
    row += [aa+'_MUT' for aa in AA]
    # trainMat += '\t'.join(['allHomologs','exclParalogs','specParalogs','orthologs','bpso','bpsh']) + '\t'
    row += ['allHomologs','exclParalogs','specParalogs','orthologs','bpso','bpsh']
    for position in range(startWS, endWS+1):
        # trainMat += ('_'+str(position)+'\t').join(['A', 'D', 'R']) + '_'+str(position)+'\t'
        row += [mut_type+'_'+str(position) for mut_type in ['A', 'D', 'R']]
    # trainMat += '\n'
    data.append(row)

    # Open the input file and convert its contents into array
    f = open(inputFile, "r")
    file_contents = f.read().split('\n')

    # Go line by line through the contents
    kinases = {}
    # for line in tqdm(open(inputFile, 'r')):
    for line in tqdm(file_contents):
        # Ignore empty line
        if line.split() == []: continue
        # Ignore comments
        if line[0] in ['#']: continue
        # Extract features
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
        # trainMat += acc + '\t' + gene + '\t' + mutation + '\t' + 'test' + '\t'
        # trainMat += str(hmmPos) + '\t' + str(hmmSS) + '\t'
        # trainMat += str(hmmScoreWT) + '\t' + str(hmmScoreMUT) + '\t' + str(hmmScoreMUT-hmmScoreWT) + '\t'
        # trainMat += str(is_phosphomimic) + '\t'
        # trainMat += '\t'.join([str(item) for item in charges_row]) + '\t'
        # trainMat += '\t'.join([str(item) for item in ptm_row]) + '\t'
        # trainMat += '\t'.join([str(item) for item in aa_row]) + '\t'
        # trainMat += '\t'.join([str(item) for item in homology_row]) + '\t'
        # trainMat += '\t'.join([str(item) for item in adr_row]) + '\n'
        
        # store in 2D array
        row = [acc, gene, mutation, 'test']
        row += [int(hmmPos), str(hmmSS)]
        row += [float(hmmScoreWT), float(hmmScoreMUT), float(hmmScoreMUT)-float(hmmScoreWT)]
        row.append(is_phosphomimic)
        row += [item for item in charges_row]
        row += [item for item in ptm_row]
        row += [item for item in aa_row]
        row += [item for item in homology_row]
        row += [item for item in adr_row]
        data.append(row)

    df = pd.DataFrame(data[1:], columns=data[0])
    
    # columns to include
    columns_to_consider = []
    for line in open(BASE_DIR+'/ML/'+'columns_to_consider.txt', 'r'):
        columns_to_consider.append(line.strip())

    # obtain the final matrix to be testes
    test_data = df.to_numpy()
    features = df.loc[:, df.columns.isin(columns_to_consider)]
    features = features.to_numpy()
    features = features[:, 3:]

    # load the prediction and scaler models
    filename = BASE_DIR+'/ML/'+'finalized_model.sav'
    clf = pickle.load(open(filename, 'rb'))
    scaler = pickle.load(open(BASE_DIR+'/ML/'+'finalized_scaler.pkl', 'rb'))
    features = scaler.transform(features)
    
    # Print the results
    results = {}
    print (''.join(['-' for i in range(50)]))
    outputText = '#Acc\tGene\tMut\tHMMpos\tPrediction\n'
    for row, predict in zip(test_data, clf.predict_proba(features)):
        outputText += row[0] +'\t'+ row[1] +'\t'+ row[2] +'\t'+ str(row[4]) +'\t'
        outputText += str(round(predict[1], 3)) + '\n'
        name = row[0] + '/' + row[2]
        results[name] = {'prediction':round(predict[1], 3), 'hmmPos':row[4]}
    if outputFile != None: open(outputFile, 'w').write(outputText)
    else: print (outputText)
    return results


if __name__ == '__main__':
    # set arguments
    parser = argparse.ArgumentParser(description='kinaseX', epilog='End of help.')
    parser.add_argument('--i', help='path to input file (mechismo-like format); see data/sample.fasta')
    parser.add_argument('--o', help='path to output file; default: print on the screen')
    args = parser.parse_args()

    # set input file
    inputFile = args.i
    if inputFile == None: inputFile = 'sample_mutations.txt'

    # set output file
    outputFile = args.o

    results_json = predict(inputFile, outputFile)