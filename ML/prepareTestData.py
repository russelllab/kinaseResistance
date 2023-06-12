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
import shutil
import gzip, sys

PTM_TYPES = ['ac', 'gl', 'm1', 'm2', 'm3', 'me', 'p', 'sm', 'ub']
MUT_TYPES = ['A', 'D', 'R']
AA = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
MODEL_NAMES = ['AIvNLD', 'LDvNAI', 'RvN', 'AIvLD']
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
    mydb = fetchData.connection(db_name='kinase_project2')
    mydb.autocommit = True
    mycursor = mydb.cursor()

    # Make header
    data = []
    row = ['Input', 'Acc','Gene','UniProtID', 'ProteinName', 'Mutation','Region', 'Dataset']
    row += ['AdjacentSites', 'alnPos', 'hmmPos','hmmSS','hmmScoreWT','hmmScoreMUT','hmmScoreDiff']
    row += ['Phosphomimic', 'Acetylmimic']
    row += ['IUPRED', 'ATPcount']
    row += ['ncontacts', 'nresidues', 'mech_intra']
    row += ['phi_psi', 'sec', 'burr', 'acc']
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
        for mut_type in ['A', 'D', 'R']:
            for aa in AA:
                row += [mut_type+'_'+aa+'_'+str(position) + '_pfam']

    for position in range(startWS, endWS+1):
        row += [mut_type+'_'+str(position) for mut_type in MUT_TYPES]
        row += [mut_type+'_'+str(position)+'_pfam' for mut_type in MUT_TYPES]
    data.append(row)

    # Open the input file and convert its contents into an array
    # f = open(inputFile, "r")
    f = fetchData.checkGZfile(inputFile)
    file_contents = f.read().split('\n')

    # Go line by line through the contents
    kinases = {}
    entries_not_found = {}
    count = 0
    for line in tqdm(file_contents):
        count += 1
        if line.split() == []: continue # Ignore empty line
        if line[0] == '#' or line.lstrip().rstrip() == '': continue # Ignore comments
        
        # Check if the line is in the correct format
        if '/' not in line:
            entries_not_found[line] = 'Incorrect format. Use per line: kinase/mutation. Eg: MAP2K1/Q56P.'
            continue

        # Check if the line has more than one /
        if line.count('/') > 1:
            entries_not_found[line] = 'Incorrect format. Use per line: kinase/mutation. Eg: MAP2K1/Q56P.'
            continue

        # Check if string before / is valid
        kinase = line.split('/')[0].rstrip().upper()
        mutation = line.split('/')[1].rstrip().upper()
        if kinase == '' or mutation == '':
            entries_not_found[line] = 'Kinase or mutation absent. Use per line: kinase/mutation. Eg: MAP2K1/Q56P.'
            continue
        
        # Check if the mutation is in the correct format
        if len(mutation) < 3:
            entries_not_found[line] = 'Incorrect mutation format. Use per line: kinase/mutation. Eg: MAP2K1/Q56P.'
            continue
        if mutation[0] not in AA:
            entries_not_found[line] = 'Incorrect wild type amino acid. '+mutation[0]+' not a valid amino acid.'
            continue
        if mutation[-1] not in AA:
            entries_not_found[line] = 'Incorrect mutatant amino acid. '+mutation[-1]+' not a valid amino acid.'
            continue
        if not mutation[1:-1].isdigit():
            entries_not_found[line] = 'Incorrect mutation format. '+mutation[1:-1]+' not a valid position. Use per line: kinase/mutation. Eg: MAP2K1/Q56P.'
            continue
        if int(mutation[1:-1]) < 1:
            entries_not_found[line] = 'Incorrect mutation format. '+mutation[1:-1]+' not a valid position. Use per line: kinase/mutation. Eg: MAP2K1/Q56P.'
            continue
        if mutation[0] == mutation[-1]:
            entries_not_found[line] = 'Incorrect mutation. Wild type and mutant amino acids are the same.'
            continue

        # When above checks are passed, run the prediction
        # Extract features i.e. name = kinase/mutation
        name = line.split()[0].rstrip().upper()
        kinase = name.split()[0].split('/')[0]
        mutation = name.split('/')[1].rstrip()

        # Retrieve acc and gene of the input kinase
        acc, gene, uniprot_id, protein_name, protein_length = fetchData.getAccGene(mycursor, kinase)

        # Check if the acc is None, which means
        # that the input kinase was not found in
        # the DB
        if acc is None:
            entries_not_found[name] = 'Protein identifier ' + kinase + ' not found. Try another identifier.'
            continue

        # Check if the mutation position is greater than the protein length
        # or if the mutation is in the first or last 2 residues
        if int(mutation[1:-1]) > protein_length:
            entries_not_found[name] = 'Position ' + str(mutation[1:-1]) +\
                                    ' is greater than the protein length'\
                                    + str(protein_length) + '.' 
            continue
        if int(mutation[1:-1]) < 3:
            entries_not_found[name] = 'Position ' + str(mutation[1:-1]) +\
                                    ' is less than 3. Given the window size of 5, ' +\
                                    'the mutation position should be greater than 2 and'+\
                                    ' less than the protein length minus 2.'
            continue
        if int(mutation[1:-1]) > protein_length - 2:
            entries_not_found[name] = 'Position ' + str(mutation[1:-1]) +\
                                    ' is greater than the protein length'\
                                    + str(protein_length) + '. Given the window size of 5, ' +\
                                    'the mutation position should be greater than 2 and'+\
                                    ' less than the protein length minus 2.'
            continue

        # Run some other checks:
        # error = 1 if the given position found in the kinase in DB
        # error = 2 if the given mutation position has the said WT AA
        error, aa_found = fetchData.checkInputPositionAA(acc, mutation, mycursor)

        # Skip if errors 1 or 2 returned
        if error == 1:
            entries_not_found[name] = 'Position ' + str(mutation[1:-1])
            entries_not_found[name] += ' not found in ' + kinase + '.'
            continue
        if error == 2:
            entries_not_found[name] = 'Did you mean ' + aa_found + mutation[1:] + '? '+\
                                    'We found a ' + aa_found + ' at position '
            entries_not_found[name] += str(mutation[1:-1]) + ' in ' + kinase
            entries_not_found[name] += ' instead of a ' + mutation[0] +\
                                    '. Please ensure that you are using the canonical '+\
                                    'isoform from UniProt (Release 2023_02).'
            continue
        
        # Make dictionary of Kinases
        # using class template
        if acc not in kinases:
            kinases[acc] = Kinase(acc, gene)
        if mutation not in kinases[acc].mutations:
            kinases[acc].mutations[mutation] = Mutation(mutation, '-', acc, 'test')

        # Get region
        region = fetchData.getRegion(mycursor, acc, mutation)
        
        # Get features
        position = int(mutation[1:-1])
        mutAA = mutation[-1]
        wtAA = mutation[0]
        hmmPos, hmmScoreWT, hmmScoreMUT, hmmSS = fetchData.getHmmPkinaseScore(mycursor, acc, wtAA, position, mutAA)
        iupred_score = fetchData.getIUPredScore(mycursor, acc, wtAA, position, mutAA)
        alnPos = fetchData.getAlnPos(mycursor, hmmPos)
        adjacentSites = fetchData.getAdjacentSites(mycursor, acc, position, 5)
        ## Even if the position is not in the HMM,
        ## we still want to show it in the output
        ptm_row = fetchData.getPTMscore(mycursor, acc, position, WS)
        aa_row = fetchData.getAAvector(wtAA, mutAA)
        homology_row = fetchData.getHomologyScores(mycursor, acc, wtAA, position, mutAA)
        if homology_row == None: continue
        mech_intra_row = fetchData.getMechIntraScores(mycursor, acc, wtAA, position, mutAA)
        if mech_intra_row == None: continue
        dssp_row = fetchData.getDSSPScores(mycursor, acc, wtAA, position, mutAA)
        if dssp_row == None: continue
        atp_count = fetchData.getATPbindingScores(mycursor, acc, position)
        is_phosphomimic = kinases[acc].mutations[mutation].checkPhosphomimic()
        is_acetylmimic = kinases[acc].mutations[mutation].checkAcetylmimic()
        charges_row = kinases[acc].mutations[mutation].findChangeInCharge()
        count_aa_change_row = fetchData.getCountAAchange(mycursor, acc, position, kinases, ws=WS)
        adr_row = fetchData.getADRvector(mycursor, acc, position, kinases, WS)
        
        # store features in 2D array
        row = [name, acc, gene, uniprot_id, protein_name, mutation, region, 'test']
        row += [adjacentSites, str(alnPos), str(hmmPos), str(hmmSS)]
        row += [float(hmmScoreWT), float(hmmScoreMUT), float(hmmScoreMUT)-float(hmmScoreWT)]
        row.append(is_phosphomimic)
        row.append(is_acetylmimic)
        row.append(iupred_score)
        row.append(atp_count)
        row += [item for item in mech_intra_row]
        row += [item for item in dssp_row]
        row += [item for item in charges_row]
        row += [item for item in ptm_row]
        row += [item for item in aa_row]
        row += [item for item in homology_row]
        row += [int(item) for item in count_aa_change_row]
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
        print (f'No data in the input file ({inputFile}) could be processed.')
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

    dic_features = {}
    dic_scalers = {}
    dic_clfs = {}
    dic_predictions = {}
    for model in MODEL_NAMES:
        dic_scalers[model] = pickle.load(open(BASE_DIR+'/ML/'+'scaler_'+model+'.pkl', 'rb'))
        dic_features[model] = dic_scalers[model].transform(features)
        dic_clfs[model] = pickle.load(open(BASE_DIR+'/ML/'+'model_'+model+'.sav', 'rb'))
        dic_predictions[model] = dic_clfs[model].predict_proba(dic_features[model])
    
    # print (dic_predictions)

    '''
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
    '''
    
    # Print the results
    terminal_width, terminal_height = shutil.get_terminal_size()
    print (''.join(['-' for i in range(terminal_width)]))
    # outputText = '#UserInput\tAcc\tGN\tUniProtID\tMutation\tHMMpos\tRegion\tPTM\tKnownADR\tNeighSites'
    # for model in MODEL_NAMES: outputText += '\t'+model
    # outputText += '\n'
    outputDict = {'UserInput': [], 'UniProtAcc': [], 'GeneName': [], 'UniProtID': [],
                  'Mutation': [], 'HMMpos': [], 'Region': [], 'PTM': [], 'KnownADR': [],
                  'NeighSites': []}
    for model in MODEL_NAMES: outputDict[model] = []
    # for row, predictAD, predictRN in zip(test_data, clfAD.predict_proba(featuresAD), clfRN.predict_proba(featuresRN)):
    for count, row in enumerate(test_data):
        ## Set prediction proba to NA if the position is not in the HMM
        user_input = row[0]
        acc = row[1]
        gene = row[2]
        uniprot_id = row[3]
        protein_name = row[4]
        mutation = row[5]
        region = row[6]
        adjacentSites = row[8]
        alnPos = row[9]
        hmmPos = row[10]
        name = user_input
        ptmType = '-'
        for ptm_type_header in [ptm_type+'_0' for ptm_type in PTM_TYPES]:
            # print (df[df['Input']==name][ptm_type_header])
            ptmTypeCell = df[df['Input']==name][ptm_type_header].values[0]
            if int(ptmTypeCell) != 0:
                ptmType = ptm_type_header.split('_')[0]
                break
        mutType = []
        mutType = fetchData.mutTypes(mycursor, acc, mutation)
        
        scores = []
        for model in MODEL_NAMES:
            if str(hmmPos) == '-':
                scores.append('NA')
            elif int(hmmPos) < 30 and int(hmmPos) > 739:
                scores.append('NA')
            else:
                scores.append(str(round(dic_predictions[model][count][1], 3)))
        # outputText += user_input +'\t'+ acc +'\t'+ gene +'\t'+ uniprot_id +'\t' + mutation +'\t'
        # outputText += str(hmmPos) +'\t'+ region + '\t' + ptmType + '\t' + ','.join(mutType) + '\t'
        # outputText += adjacentSites + '\t'
        # outputText += '\t'.join(scores) + '\n'
        outputDict['UserInput'].append(user_input)
        outputDict['UniProtAcc'].append(acc)
        outputDict['GeneName'].append(gene)
        outputDict['UniProtID'].append(uniprot_id)
        outputDict['Mutation'].append(mutation)
        outputDict['HMMpos'].append(hmmPos)
        outputDict['Region'].append(region)
        outputDict['PTM'].append(ptmType)
        outputDict['KnownADR'].append(','.join(mutType))
        outputDict['NeighSites'].append(adjacentSites)
        for model, score in zip(MODEL_NAMES, scores):
            outputDict[model].append(score)
        if len(mutType) == 0: mutType = '-'
        else: mutType = ''.join(mutType)
        results['predictions'][name] = {
                        'acc':acc,
                        'gene':gene,
                        'uniprot_id': uniprot_id,
                        'protein_name': protein_name,
                        'mutation':mutation,
                        'adjacentSites': adjacentSites,
                        # 'predAD':prediction_probAD,
                        # 'predRN':prediction_probRN,
                        'hmmPos':hmmPos,
                        'alnPos':alnPos,
                        'region':region,
                        'ptmType':ptmType,
                        'mutType':mutType,
                        }
        for model in MODEL_NAMES:
            if str(hmmPos) == '-':
                results['predictions'][name][model] = 'NA\t'
            elif int(hmmPos) < 30 and int(hmmPos) > 739:
                results['predictions'][name][model] = 'NA\t'
            else:
                results['predictions'][name][model] = str(round(dic_predictions[model][count][1], 3)) + '\t'
        
            # if 'A84' in mutation and model=='AIvNLD':
                # print (dic_features[model][count])
                # print (len(dic_features[model][count]))
                # for f, hh in zip(columns_to_consider, dic_features[model][count]):
                #     print (f, hh)
    
    outputDF = pd.DataFrame(outputDict)
    if outputFile != None: gzip.open(outputFile+'.gz', 'wt').write(outputDF.to_string(index=False))
    # else: print (outputText)
    else: print (outputDF.to_string(index=False))
    return results
    # yield results

# Run this script from command line
if __name__ == '__main__':
    # set arguments
    parser = argparse.ArgumentParser(description='Activaark', epilog='End of help.')
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
