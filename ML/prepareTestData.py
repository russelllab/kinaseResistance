#!/usr/bin/env python3.10
# coding: utf-8

"""
This script is used to prepare the test data for the kinase activity prediction
and make the predictions. The input file should be a text file with the following
format:
    MAP2K1/Q56P
    PIM1/T23I
    P11309/S97N
    P11309/Q127E
"""

import fetchData
from tqdm import tqdm
import pickle
import pandas as pd
from cls import Kinase, Mutation
import argparse

PTM_TYPES = ['ac', 'gl', 'm1', 'm2', 'm3', 'me', 'p', 'sm', 'ub']
MUT_TYPES = ['A', 'D', 'R']
AA = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
WS = 5

def predict(inputFile, outputFile = None, BASE_DIR = '../') -> dict:
    """
    Function to predict the effect of mutations on kinase activity

    Parameters:
    inputFile (str): Path to the input file
    outputFile (str): Path to the output file
    BASE_DIR (str): Path to the base directory

    Returns:
    results (dict): Dictionary containing the results

    Contents of the results dictionary:
    entries_not_found (dict): Dictionary containing the entries that
    were not found in the database and the reason why
    predictions (dict): Dictionary containing the predictions for each entry
    and some other information
    """

    # Create dict results that will return the output
    results = {
            'entries_not_found': {},
            'predictions': {}
            }

    # Connect to DB
    mydb = fetchData.connection()
    mydb.autocommit = True
    mycursor = mydb.cursor()

    # Make header
    data = []
    row = ['Input', 'Acc','Gene','Mutation','Dataset']
    row += ['hmmPos','hmmSS','hmmScoreWT','hmmScoreMUT','hmmScoreDiff']
    row += ['Phosphomimic']
    row += ['ChargesWT','ChargesMUT','ChargesDiff']
    startWS = int((WS-1)/2) * -1
    endWS = int((WS-1)/2)
    for position in range(startWS, endWS+1):
        row += [ptm_type+'_'+str(position) for ptm_type in PTM_TYPES]
        row += [ptm_type+'_'+str(position)+'_pfam' for ptm_type in PTM_TYPES]
    row += [aa+'_WT' for aa in AA]
    row += [aa+'_MUT' for aa in AA]
    row += ['allHomologs','exclParalogs','specParalogs','orthologs','bpso','bpsh']
    for position in range(startWS, endWS+1):
        row += [mut_type+'_'+str(position) for mut_type in MUT_TYPES]
        row += [mut_type+'_'+str(position)+'_pfam' for mut_type in MUT_TYPES]
    data.append(row)

    # Open the input file and convert its contents into an array
    f = open(inputFile, "r")
    file_contents = f.read().split('\n')

    # Go line by line through the contents
    kinases = {}
    entries_not_found = {}
    count = 0
    for line in tqdm(file_contents):
        count += 1
        if line.split() == []: continue # Ignore empty line
        if line[0] == '#' or line.lstrip().rstrip() == '': continue # Ignore comments
        
        # Extract features i.e. name = kinase/mutation
        name = line.split()[0].rstrip().upper()
        kinase = name.split()[0].split('/')[0]
        mutation = name.split('/')[1].rstrip()

        # Retrieve acc and gene of the input kinase
        acc, gene, uniprot_id = fetchData.getAccGene(mycursor, kinase)

        # Check if the acc is None, which means
        # that the input kinase was not found in
        # the DB
        if acc is None:
            entries_not_found[name] = 'Protein identifier ' + kinase + ' not found'
            continue

        # Run some other checks:
        # error = 1 if the given position found in the kinase in DB
        # error = 2 if the given mutation position has the said WT AA
        error, aa_found = fetchData.checkInputPositionAA(acc, mutation, mycursor)

        # Skip if errors 1 or 2 returned
        if error == 1:
            entries_not_found[name] = 'Position ' + str(mutation[1:-1])
            entries_not_found[name] += ' not found in ' + kinase
            continue
        if error == 2:
            entries_not_found[name] = 'Found ' + aa_found + ' at position '
            entries_not_found[name] += str(mutation[1:-1]) + ' in ' + kinase
            entries_not_found[name] += ' instead of ' + mutation[0]
            continue
        
        # Make dictionary of Kinases
        # using class template
        if acc not in kinases:
            kinases[acc] = Kinase(acc, gene)
        if mutation not in kinases[acc].mutations:
            kinases[acc].mutations[mutation] = Mutation(mutation, '-', acc, 'test')
        
        # Get features
        position = int(mutation[1:-1])
        mutAA = mutation[-1]
        wtAA = mutation[0]
        hmmPos, hmmScoreWT, hmmScoreMUT, hmmSS = fetchData.getHmmPkinaseScore(mycursor, acc, wtAA, position, mutAA)
        ## Even if the position is not in the HMM,
        ## we still want to show it in the output
        ptm_row = fetchData.getPTMscore(mycursor, acc, position, WS)
        aa_row = fetchData.getAAvector(wtAA, mutAA)
        homology_row = fetchData.getHomologyScores(mycursor, acc, wtAA, position, mutAA)
        if homology_row == None: continue
        is_phosphomimic = kinases[acc].mutations[mutation].checkPhosphomimic()
        charges_row = kinases[acc].mutations[mutation].findChangeInCharge()
        adr_row = fetchData.getADRvector(mycursor, acc, position, kinases, WS)
        
        # store features in 2D array
        row = [name, acc, gene, mutation, 'test']
        row += [str(hmmPos), str(hmmSS)]
        row += [float(hmmScoreWT), float(hmmScoreMUT), float(hmmScoreMUT)-float(hmmScoreWT)]
        row.append(is_phosphomimic)
        row += [item for item in charges_row]
        row += [item for item in ptm_row]
        row += [item for item in aa_row]
        row += [item for item in homology_row]
        row += [item for item in adr_row]
        # print (row)
        data.append(row)
        # yield count/float(len(file_contents)) * 100

    # save entries not found in the results dic
    results ['entries_not_found'] = entries_not_found
    
    # convert the 2D array into dataframe
    df = pd.DataFrame(data[1:], columns=data[0])
    # print (df[['Input', 'A_0', 'D_0', 'R_0']])

    # if no values in data (besides the header)
    # then just end it here and return the dic results
    if len(data) == 1:
        print ('No data found in the input file.')
        return results
        # yield results
    # else go ahead

    # columns to consider
    columns_to_consider = []
    for line in open(BASE_DIR+'/ML/'+'columns_to_consider.txt', 'r'):
        columns_to_consider.append(line.strip())

    # obtain the final matrix to be testes
    test_data = df.to_numpy()
    features = df.loc[:, df.columns.isin(columns_to_consider)]
    features = features.to_numpy()
    features = features[:, 3:]

    # load the AD prediction and scaler models
    # and scale the features matrix
    filenameAD = BASE_DIR+'/ML/'+'finalized_model_AD.sav'
    clfAD = pickle.load(open(filenameAD, 'rb'))
    scalerAD = pickle.load(open(BASE_DIR+'/ML/'+'finalized_scaler_AD.pkl', 'rb'))
    featuresAD = scalerAD.transform(features)

    # load the RN prediction and scaler models
    # and scale the features matrix
    filenameRN = BASE_DIR+'/ML/'+'finalized_model_RN.sav'
    clfRN = pickle.load(open(filenameRN, 'rb'))
    scalerRN = pickle.load(open(BASE_DIR+'/ML/'+'finalized_scaler_RN.pkl', 'rb'))
    featuresRN = scalerRN.transform(features)
    
    # Print the results
    print (''.join(['-' for i in range(50)]))
    outputText = '# Input\tAcc\tGene\tMut\tHMMpos\tPredAD\tPredRN\n'
    for row, predictAD, predictRN in zip(test_data, clfAD.predict_proba(featuresAD), clfRN.predict_proba(featuresRN)):
        ## Set prediction proba to NA if the position is not in the HMM
        acc = row[1]
        gene = row[2]
        mutation = row[3]
        if str(row[5]) != '-':
            prediction_probAD = round(predictAD[1], 3)
            prediction_probRN = round(predictRN[1], 3)
        else:
            prediction_probAD = 'NA'
            prediction_probRN = 'NA'
        outputText += row[0] +'\t'+ row[1] +'\t'+ row[2] +'\t'+ row[3] +'\t'+ str(row[5]) +'\t'
        outputText += str(prediction_probAD) + '\t' +str(prediction_probRN) + '\n'
        name = row[0]
        ptmType = '-'
        for ptm_type_header in [ptm_type+'_0' for ptm_type in PTM_TYPES]:
            # print (df[df['Input']==name][ptm_type_header])
            ptmTypeCell = df[df['Input']==name][ptm_type_header].values[0]
            if int(ptmTypeCell) != 0:
                ptmType = ptm_type_header.split('_')[0]
                break
        mutType = []
        for mut_type_header in [mut_type+'_0' for mut_type in MUT_TYPES]:
            mutTypeCell = df[df['Input']==name][mut_type_header].values[0]
            if int(mutTypeCell) != 0:
                # mutType = mut_type_header.split('_')[0]
                mutType.append(mut_type_header.split('_')[0])
                # break
        if len(mutType) == 0: mutType = '-'
        else: mutType = ''.join(mutType)
        results['predictions'][name] = {
                        'acc':acc,
                        'gene':gene,
                        'mutation':mutation,
                        'predAD':prediction_probAD,
                        'predRN':prediction_probRN,
                        'hmmPos':row[5],
                        'ptmType':ptmType,
                        'mutType':mutType,
                        }
    # print (results)
    if outputFile != None: open(outputFile, 'w').write(outputText)
    else: print (outputText)
    return results
    # yield results

# Run this script from command line
if __name__ == '__main__':
    # set arguments
    parser = argparse.ArgumentParser(description='kinaseX', epilog='End of help.')
    parser.add_argument('i', help='path to input file (mechismo-like format); see sample_mutations.txt')
    parser.add_argument('--o', help='path to output file; default: print on the screen')
    args = parser.parse_args()

    # set input file to default if not provided
    inputFile = args.i
    if inputFile == None: inputFile = 'sample_mutations.txt'

    # extract output file if provided
    outputFile = args.o

    results_json = predict(inputFile, outputFile)
    if len(results_json['entries_not_found']) > 0:
        if len(results_json['entries_not_found']) == 1: print ('The following input was ignored:')
        else: print ('The following inputs were ignored:')
        for name in results_json['entries_not_found']:
            print (name+':', results_json['entries_not_found'][name])