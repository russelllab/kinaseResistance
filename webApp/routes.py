#!/usr/bin/env python3
'''
This script contains functions called by JS and API endpoints
Routes are defined as AJAX<name> eg: AJAXvariants
Functions of routes are defined as get_<name> eg: get_variants
'''
from multiprocessing import Value
from flask import Flask, render_template, request, jsonify, redirect, url_for
import gzip, json, random, string, os, requests, urllib
from urllib.parse import urlparse
from urllib.parse import parse_qs
import re, sys, math, ast
import time
import json
import zlib
from urllib.parse import urlparse, parse_qs, urlencode
import requests
from requests.adapters import HTTPAdapter, Retry
import socket
import psycopg2

# Set up the base_url
if socket.gethostname() == 'pevolution2.bioquant.uni-heidelberg.de':
	BASE_URL = 'http://activark.russelllab.org/'
	# BASE_DIR = '/net/home.isilon/ds-russell/kinaseResistance/'
	BASE_DIR = '/var/www/flask_apps/kinaseResistance/'
elif socket.gethostname() == 'gd':
	BASE_URL = 'http://127.0.0.1:5000/'
	BASE_DIR = '/home/pevolution/projects/kinaseResistance/'
else:
	BASE_URL = 'http://127.0.0.1:5000/'
	BASE_DIR = '../'
	BASE_DIR = '/home/gurdeep/projects/kinaseResistance/'

sys.path.insert(1, BASE_DIR+'/ML/')
import prepareTestData
import fetchData

create_svg_path = '/Create_SVG/Enhancements_May2023/27June/'
sys.path.insert(1, BASE_DIR+create_svg_path)
import create_svg_20230706_kinases_GS as create_svg
conservation_dic_path = BASE_DIR+create_svg_path+'GenerelleKonservierung_Jun-28-2023.txt'
identity_dic_path = BASE_DIR+create_svg_path+'SeqIdentity_Matrix_Jun-28-2023.txt'

def connection(database='kinase_project2'):
    '''Function to connect to postgresql database'''
    mydb = psycopg2.connect(
                            database = database,
                            user = "gurdeep",
                            password = "hellokitty",
                            host = "localhost",
                            port = "5432")
    
    mydb.autocommit = True
    mycursor = mydb.cursor()
    return mycursor

def extract_pubmed_ids(string):
	# https://pubmed.ncbi.nlm.nih.gov/11781872/
	pattern = r"PubMed:\s*(\d+)"
	matches = re.findall(pattern, string)
	return [str(match) for match in matches]

def makeWindowText(position_of_interest, current_position):
	window = int(position_of_interest) - int(current_position)
	if window > 0: window = '+'+str(window)
	else: window = str(window)
	return '('+window+')'

def makeText(acc, gene, mutation, interested_kinase_pfampos, mycursor):
	'''
	Make text for prediction
	'''
	ws = 5
	if ws > 0: ws -= 1
	ws = int(ws/2)
	dic_mutations = {'A': 'Activating', 'D': 'Deactivating', 'R': 'Resistance'}
	dic_mutations = {'activating': 'Activating', 'increase': 'Activating',\
		  			'loss': 'Deactivating', 'decrease': 'Deactivating',\
					'resistance': 'Resistance'}
	dic_ptms = {'p': 'Phosphorylation-site', 'ub': 'Ubiquitination-site',
	     		'ac': 'Acetylation-site', 'me': 'Methylation-site',\
				'gl': 'Glycosylation-site', 'sm': 'Sumoylation-site',\
				'm1': 'Myristoylation-site', 'm2': 'Palmitoylation-site',\
				'm3': 'Myristoylation-site'}
	text = ''
	mutation_position = int(mutation[1:-1])
	mycursor.execute("SELECT alnpos FROM hmm \
						WHERE pfampos = %s", (str(interested_kinase_pfampos), ))
	interested_kinase_alnpos = mycursor.fetchone()[0]
	# print (interested_kinase_alnpos, interested_kinase_pfampos, gene)
	data = []
	# If the pfampos does not exist
	if interested_kinase_pfampos=='-':
		return data, text
	for position in range(mutation_position-ws, mutation_position+ws+1):
	# 	# do it for the acc of interest
	# 	mycursor.execute("SELECT pfampos FROM positions \
	# 		WHERE acc = %s and uniprotpos = %s", (acc, str(position)))
	# 	pfampos = mycursor.fetchone()[0]
	# 	mycursor.execute("SELECT uniprotaa, ptmtype FROM ptms \
	# 					WHERE acc = %s and uniprotpos = %s", (acc, str(position)))
	# 	hits = mycursor.fetchall()
	# 	for entry in hits:
	# 		text += "<b>"+str(entry[0]) + str(position)+ "</b>" + ' is a known ' + dic_ptms[entry[1]] + ' site. '
	# 		text += '<a href=\"http://www.phosphosite.org/uniprotAccAction?id='+ acc +'\" target=\"_blank\">PhosphoSitePlus <i class="bi bi-box-arrow-in-up-right"></i></a>'
	# 		text += '<br>' # add a full stop at the end of the sentence
	# 		row = []
	# 		row.append('-')
	# 		row.append(gene)
	# 		row.append(acc)
	# 		row.append(entry[0]+str(position))
	# 		row.append(pfampos)
	# 		row.append(entry[1])
	# 		row.append(dic_ptms[entry[1]])
	# 		if gene == 'BRAF': print ('gene', row)
	# 		row.append('<a href=\"http://www.phosphosite.org/uniprotAccAction?id='+ acc +'\" target=\"_blank\">PhosphoSitePlus <i class="bi bi-box-arrow-in-up-right"></i></a>')
	# 		data.append(row)
		
		# do it for all the accs at the position of interest
		mycursor.execute("SELECT alnpos, pfampos FROM positions \
			WHERE acc = %s and uniprotpos = %s", (acc, str(position)))
		alnpos, pfampos = mycursor.fetchone()
		# print(pfampos, interested_kinase_pfampos, alnpos, interested_kinase_alnpos)
		if pfampos == '-' and alnpos == '-': continue
		mycursor.execute("SELECT uniprotaa, uniprotpos, ptmtype, acc, gene FROM ptms \
						WHERE pfampos = %s", (str(pfampos), ))
		hits = mycursor.fetchall()
		for entry in hits:
			uniprotaa = entry[0]
			uniprotpos = entry[1]
			ptmtype = entry[2]
			ref_acc = entry[3]
			# if ref_acc == acc: continue
			ref_gene = entry[4]
			text += "<b>"+ref_gene+'/'+str(uniprotaa)+str(uniprotpos)+"</b>" +' is a known ' + dic_ptms[ptmtype] + ' site. '
			text += '<a href=\"http://www.phosphosite.org/uniprotAccAction?id='+ ref_acc +'\" target=\"_blank\">PhosphoSitePlus <i class="bi bi-box-arrow-in-up-right"></i></a>'
			text += '<br>' # add a full stop at the end of the sentence
			row = []
			# row.append('-')
			row.append(ref_gene)
			# row.append(ref_acc)
			row.append('<a href=\"https://www.uniprot.org/uniprot/'+ref_acc+'\" target=\"_blank\">'+ref_acc+'<i class="bi bi-box-arrow-in-up-right"></i></a>')
			row.append(uniprotaa)
			row.append(uniprotpos)
			row.append('-') # MUT aa = blank for PTMs
			# row.append(str(pfampos) + makeWindowText(pfampos, interested_kinase_pfampos))
			# row.append(str(alnpos) + makeWindowText(alnpos, interested_kinase_alnpos))
			row.append(str(position-mutation_position))
			row.append(dic_ptms[ptmtype])
			row.append('-') # Description = blank for PTMs
			# if ref_gene == 'BRAF': print ('ref_gene', row)
			row.append('<a href=\"http://www.phosphosite.org/uniprotAccAction?id='+ ref_acc +'\" target=\"_blank\">PhosphoSitePlus <i class="bi bi-box-arrow-in-up-right"></i></a>')
			data.append(row)

	# print (data)
	
	for position in range(mutation_position-ws, mutation_position+ws+1):
		# do it for acc of interest
		# mycursor.execute("SELECT pfampos FROM positions \
		# 	WHERE acc = %s and uniprotpos = %s", (acc, str(position)))
		# pfampos = mycursor.fetchone()[0]
		# mycursor.execute("SELECT mutation, mut_type, info FROM mutations \
		# 				WHERE acc = %s and wtpos = %s", (acc, str(position)))
		# hits = mycursor.fetchall()
		# for entry in hits:
		# 	mutation = entry[0]
		# 	mut_type = entry[1]
		# 	info = entry[2]
		# 	text += "<b>" + str(mutation) + "</b>" + ' is a known '+dic_mutations[mut_type]+' mutation.'
		# 	row = []
		# 	row.append('-')
		# 	row.append(gene)
		# 	row.append(acc)
		# 	row.append(mutation)
		# 	row.append(pfampos)
		# 	row.append(mut_type)
		# 	row.append(dic_mutations[mut_type])
		# 	if mut_type != 'R':
		# 		text += ' <u>Description</u>: ' + info.split('"""')[0]
		# 		pubmed_ids = extract_pubmed_ids(info.replace('"', '')) # remove double quotes
		# 		pubmed_ids_text = []
		# 		for pubmed_id in pubmed_ids:
		# 			pubmed_ids_text.append('<a href=\"https://pubmed.ncbi.nlm.nih.gov/' + str(pubmed_id) + '\" target=\"_blank\">' + str(pubmed_id) + '<i class="bi bi-box-arrow-in-up-right"></i></a>')
		# 		if pubmed_ids_text != []:
		# 			text += ' <u>PubMed</u>: ' + '; '.join(pubmed_ids_text)
		# 			row.append('<u>PubMed</u>: ' + '; '.join(pubmed_ids_text))
		# 		else:
		# 			row.append('-')
		# 		text += '.<br>' # add a full stop at the end of the sentence
		# 	else:
		# 		cosmic = ' ' + '<a href=\"https://cancer.sanger.ac.uk/cosmic/gene/analysis?ln='
		# 		cosmic += gene +'#drug-resistance\" target=\"_blank\">COSMIC <i class="bi bi-box-arrow-in-up-right"></i></a>'
		# 		text += cosmic
		# 		row.append(cosmic)
		# 	data.append(row)
		
		# do it for all the accs at the position of interest
		mycursor.execute("SELECT alnpos, pfampos FROM positions \
			WHERE acc = %s and uniprotpos = %s", (acc, str(position)))
		alnpos, pfampos = mycursor.fetchone()
		# print(pfampos)
		if pfampos == '-' and alnpos == '-': continue
		mycursor.execute("SELECT mutation, wtaa, wtpos, mut_type, acc, gene, info, pubmed FROM mutations \
						WHERE pfampos = %s", (str(pfampos), ))
		hits = mycursor.fetchall()
		for entry in hits:
			ref_mutation = entry[0]
			uniprotaa = entry[1]
			uniprotpos = entry[2]
			mut_type = entry[3]
			if mut_type == 'neutral': continue
			ref_acc = entry[4]
			# skip the acc of interest since it has been already done
			# if ref_acc == acc: continue
			ref_gene = entry[5]
			info = entry[6]
			pubmedIDs = entry[7]
			text += "<b>" + ref_gene+'/'+str(ref_mutation) + "</b>" + ' is a known '+dic_mutations[mut_type]+' mutation.'
			row = []
			row.append(ref_gene)
			# row.append(ref_acc)
			row.append('<a href=\"https://www.uniprot.org/uniprot/'+ref_acc+'\" target=\"_blank\">'+ref_acc+'<i class="bi bi-box-arrow-in-up-right"></i></a>')
			row.append(ref_mutation[0])
			row.append(ref_mutation[1:-1])
			row.append(ref_mutation[-1])
			# row.append(str(pfampos) + makeWindowText(pfampos, interested_kinase_pfampos))
			# row.append(str(alnpos) + makeWindowText(alnpos, interested_kinase_alnpos))
			row.append(str(position-mutation_position))
			row.append(dic_mutations[mut_type])
			# row.append(info.split('"""')[0] if '"' in info else '-')
			row.append(info)
			if mut_type != 'resistance':
				text += ' <u>Description</u>: ' + info.split('"""')[0]
				# pubmed_ids = extract_pubmed_ids(info.replace('"', '')) # remove double quotes
				pubmed_ids = pubmedIDs.split(',')
				pubmed_ids_text = '+or+'.join(pubmed_ids)
				# for pubmed_id in pubmed_ids:
					# pubmed_ids_text.append('<a href=\"https://pubmed.ncbi.nlm.nih.gov/' + str(pubmed_id) + '\" target=\"_blank\">' + str(pubmed_id) + '<i class="bi bi-box-arrow-in-up-right"></i></a>')
				if pubmed_ids != []:
					# text += ' <u>PubMed</u>: ' + '; '.join(pubmed_ids_text)
					text += '<a href=\"https://pubmed.ncbi.nlm.nih.gov/?term=' + str(pubmed_ids_text) + '\" target=\"_blank\">PubMed<i class="bi bi-box-arrow-in-up-right"></i></a>'
					# row.append('<u>PubMed</u>: ' + '; '.join(pubmed_ids_text))
					row.append('<a href=\"https://pubmed.ncbi.nlm.nih.gov/?term=' + str(pubmed_ids_text) + '\" target=\"_blank\">PubMed<i class="bi bi-box-arrow-in-up-right"></i></a>')
				else:
					row.append('-')
				text += '.<br>' # add a full stop at the end of the sentence
			else:
				cosmic = ' ' + '<a href=\"https://cancer.sanger.ac.uk/cosmic/gene/analysis?ln='
				cosmic += ref_gene +'#drug-resistance\" target=\"_blank\">COSMIC <i class="bi bi-box-arrow-in-up-right"></i></a>'
				text += cosmic
				row.append(cosmic)
			data.append(row)
	# print (text)
	return data, text

def makeUniqID():
	'''
	Generates a unique ID for the job
	'''
	# A unique ID is created before the pocess begins
	# Each submitted job is assigned a unique ID
	uniqID = ''.join(random.choices(string.ascii_uppercase + string.digits, k=5))
	return uniqID

def resetDic(dic, alignment):
	'''
	dic[acc][ptm_type_name] = [str(position), wtAA+str(position), text]
	dic[acc][mut_type_name] = [str(position), wtAA+str(position)+mutAA, text]
	'''
	new_dic = {}
	for line in open(alignment, 'r'):
		if line.startswith('#') or line.startswith('\n'): continue
		name = line.split()[0]
		acc = line.split()[0].split('|')[1]
		if acc not in dic: continue
		# print (line)
		start, end = line.split()[0].split('|')[2].split('-')
		for type_name in dic[acc]:
			for row in dic[acc][type_name]:
				position, site, text = row
				if int(position) >= int(start) and int(position) <= int(end):
					if name not in new_dic: new_dic[name] = {}
					if type_name not in new_dic[name]: new_dic[name][type_name] = []
					new_dic[name][type_name].append(row)
	# geeky_file = open('sample_dic_mutation_info.txt', 'wt')
	# geeky_file.write(str(new_dic))
	# print (new_dic)
	return new_dic

def final_verdictNDA(AIvLD, A, D, N, RvN):
	if AIvLD == 'NA\t':
		verdict = 'NA'
		return verdict
	'''
	a = float(A)*0.6 + float(AIvLD)*0.4
	d = float(D)*0.6 + (1-float(AIvLD))*0.4
	n = float(N)*0.6
	if n > 0.5:
		verdict = '-'
		return verdict
	if a >= d:
		prob = a
		verdict = 'Activating'
	else:
		prob = d
		verdict = 'Deactivating'
	if 0 <= prob < 0.33:
		verdict += ' (Low)'
	elif 0.33 <= prob < 0.66:
		verdict += ' (Medium)'
	else:
		verdict += ' (High)'
	return verdict
	'''
	verdictA = '-'
	if float(AIvLD) >= 0.7 and float(A) > float(D) and float(A) > float(N):
		verdictA = 'Activating (High)'
	elif float(AIvLD) >= 0.5 and float(A) > float(D) and float(A) > float(N):
		verdictA = 'Activating (Medium)'
	elif float(AIvLD) >= 0.5 and float(N) > float(A) and float(A) > float(D):
		verdictA = 'Activating (Low)'
	
	verdictD = '-'
	if float(AIvLD) <= 0.3 and float(D) > float(A) and float(D) > float(N):
		verdictD = 'Deactivating (High)'
	elif float(AIvLD) < 0.5 and float(D) > float(A) and float(D) > float(N):
		verdictD = 'Deactivating (Medium)'
	elif float(AIvLD) < 0.5 and float(N) > float(D) and float(D) > float(A):
		verdictD = 'Deactivating (Low)'
	
	if verdictA == '-' and verdictD == '-':
		verdict = 'Uncertain'
	elif verdictA != '-' and verdictD == '-':
		verdict = verdictA
	elif verdictA == '-' and verdictD != '-':
		verdict = verdictD
	else:
		if float(AIvLD) >= 0.6:
			verdict = verdictA
		elif float(AIvLD) <= 0.4:
			verdict = verdictD
		else:
			verdict = 'Uncertain'
	
	return verdict
	
def final_verdictR(AIvLD, A, D, RvN):
	if RvN == 'NA\t':
		verdict = '-'
		return verdict
	if float(RvN) > 0.5 and float(AIvLD) > 0.3:
		verdict = 'Resistance'
	else:
		verdict = '-'
	return verdict


def makeOutputJson(uniqID, results, mycursor) -> dict:
	output = []
	for name in results['predictions']:
		# num += 1
		# name = acc + '/' + mutation
		# print (name)
		kinase, mutation = name.split('/')
		dic = {}
		dic['name'] = name
		dic['mutation'] = results['predictions'][name]['mutation']
		# https://www.uniprot.org/uniprotkb/P05067/entry
		# dic['acc'] = results['predictions'][name]['acc']
		acc = results['predictions'][name]['acc']
		dic['acc'] = acc
		# dic['acc'] = '<a href=\"https://www.uniprot.org/uniprot/'+acc+'\" target=\"_blank\">'+acc+'<i class="bi bi-box-arrow-in-up-right"></i></a>'
		dic['view'] = '<a href=\"'
		dic['view'] += BASE_URL
		dic['view'] += 'result?uniqID='+uniqID+'&kinase=' + kinase + '&mutation=' + mutation
		# dic['view'] += '\" target=\"_blank\">View</a>'
		dic['view'] += '\" target=\"_blank\"><i class="bi bi-box-arrow-in-up-right"></i></a>'
		# dic['view'] += '\" target=\"_blank\"></a>'
		# dic['view'] += '<a href="#"><i class="bi bi-box-arrow-in-up-right"></i></a>'
		dic['gene'] = results['predictions'][name]['gene']
		dic['uniprot_id'] = results['predictions'][name]['uniprot_id']
		dic['protein_name'] = results['predictions'][name]['protein_name']
		dic['region'] = results['predictions'][name]['region']
		dic['ptmType'] = results['predictions'][name]['ptmType']
		dic['mutType'] = results['predictions'][name]['mutType']
	
		dic['AIvNLD'] = results['predictions'][name]['AIvNLD']
		dic['LDvNAI'] = results['predictions'][name]['LDvNAI']
		dic['RvN'] = results['predictions'][name]['RvN']
		dic['AIvLD'] = results['predictions'][name]['AIvLD']
		# print (dic['AIvLD'].split(), name)
		if dic['AIvLD'].lstrip().rstrip() == 'NA':
			# print (name)
			dic['view'] = '-'
		dic['AIvN'] = results['predictions'][name]['AIvN']
		dic['LDvN'] = results['predictions'][name]['LDvN']
		dic['RvN'] = results['predictions'][name]['RvN']

		# NDA predictions
		dic['A'] = results['predictions'][name]['A']
		dic['D'] = results['predictions'][name]['D']
		dic['N'] = results['predictions'][name]['N']

		# Verdict
		dic['verdictNDA'] = final_verdictNDA(
										dic['AIvLD'],
				 						dic['A'],
										dic['D'],
										dic['N'],
										dic['RvN']
										)
		
		dic['verdictR'] = final_verdictR(
										dic['AIvLD'],
				 						dic['A'],
										dic['D'],
										dic['RvN']
										)
		
		dic['hmmPos'] = results['predictions'][name]['hmmPos']
		dic['alnPos'] = results['predictions'][name]['alnPos']
		adjacentSites = results['predictions'][name]['adjacentSites']
		centerAdjacentSites = (len(adjacentSites)-1)/2
		adjacentSitesBold = ''
		for num, char in enumerate(adjacentSites):
			if num == centerAdjacentSites:
				adjacentSitesBold += '<b><u>' + char + '</u></b>'
			else:
				adjacentSitesBold += char
		# dic['adjacentSites'] = results['predictions'][name]['adjacentSites']
		dic['adjacentSites'] = adjacentSitesBold
		_, dic['text'] = makeText(dic['acc'], dic['gene'], dic['mutation'], dic['hmmPos'], mycursor)
		output.append(dic)
		# yield str(num) + '\n'
	
	output = {'data': output}
	return output

def runPrediction(uniqID, inputMuts):
	if os.path.isfile(BASE_DIR + '/webApp/static/predictor/output/'+uniqID) is False:
		os.system('mkdir '+ BASE_DIR +'/webApp/static/predictor/output/'+uniqID)
	with open(BASE_DIR+'/webApp/static/predictor/output/'+uniqID+'/input.txt', 'w') as f:
		f.write(inputMuts)

	# progress = prepareTestData.predict(BASE_DIR+'/webApp/static/predictor/output/'+uniqID+'/input.txt', \
	# 			BASE_DIR = BASE_DIR)
	# for p in progress:
	# 	if isinstance(p, float):
	# 		print(f"Task is {p}% complete")
	# 	else:
	# 		print(p.keys())
	# 		return p
		# yield p
	results = prepareTestData.predict(10, BASE_DIR+'/webApp/static/predictor/output/'+uniqID+'/input.txt', \
				BASE_DIR = BASE_DIR)
	return results

def configureRoutes(app):
	# print (os.getcwd())
	alignmentsPath = BASE_DIR + '/analysis/alignments/data/HUMAN/orthologs_only/'
	mutationsPath = BASE_DIR + '/analysis/mutations/data/VLatest/'
	featuresPath = BASE_DIR + '/analysis/features/data/VLatest/'
	interactorsPath = BASE_DIR + '/analysis/interactors/data/VLatest/'
	speciesTaxonomyPath = BASE_DIR + '/analysis/tax_lineage_table/'
	idmappingPath = BASE_DIR + '/webApp/static/data/'

	@app.route('/', methods=['GET', 'POST'])
	def home():
		''''
		Generates a unique ID for the job and renders the home page
		'''
		# A unique ID is created before the pocess begins
		# Each submitted job is assigned a unique ID
		# uniqID = ''.join(random.choices(string.ascii_uppercase + string.digits, k=5))
		# return render_template('maintenance.html')
		uniqID = makeUniqID()
		# return 'Home'
		# return render_template('home.html', uniqID=uniqID)
		return render_template('progress2.html', uniqID=uniqID)
	
	@app.route('/alignment', methods=['GET', 'POST'])
	def alignment():
		''''
		Display the alignment page
		'''
		return render_template('alignment.html')
	
	@app.route('/about', methods=['GET', 'POST'])
	def about():
		''''
		Display the about page
		'''
		return render_template('about.html')
	
	@app.route('/help', methods=['GET', 'POST'])
	def help():
		''''
		Display the help page
		'''
		return render_template('help.html')
	
	@app.route('/datasets', methods=['GET', 'POST'])
	def datasets():
		''''
		Display the datasets page
		'''
		return render_template('datasets.html')
	
	@app.route('/examples', methods=['GET', 'POST'])
	def examples():
		''''
		Display the examples page
		'''
		return render_template('examples.html')

	@app.route('/output/<string:uniqID>', methods=['GET', 'POST'])
	def output(uniqID: str):
		'''
		A route executed when the user clicks the submit button
		on the home page. It takes uniqID as the parameter. Creates
		a summary data table with mutation and related information.
		'''
		if request.method == 'POST':
			inputMuts = request.form['inputMut']
		if request.method == 'GET':
			# inputMuts = 'BRAF/V600E'
			# f = open(BASE_DIR+'/tests/sample_mutations3.txt', "r")
			f = open(BASE_DIR+'/webApp/static/predictor/output/'+uniqID+'/input.txt', "r")
			inputMuts = f.read()

		mycursor = connection()
		#sys.exit()
		# inputMuts = request.form['inputMut']
		results = runPrediction(uniqID, inputMuts)
		# print (results)
		if isinstance(results, float):
			return app.response_class(results, mimetype='text/plain')

		# if len(results['predictions']) == 0:
		# 	return render_template('error1.html',
		# 						# flaggedInput=json.dumps(kinase+'/'+mutation)
		# 						)
		
		ignored = []
		entries_not_found = results['entries_not_found']
		print (entries_not_found)
		for name in results['entries_not_found']:
			# kinase, mutation = name.split('/')
			dic = {}
			dic['name'] = name
			dic['reason'] = results['entries_not_found'][name]
			ignored.append(dic)
		ignored = {'data': ignored}

		# error_html_text = ''
		# for name in entries_not_found:
		# 	error_html_text += name + ' ' + entries_not_found[name] + '<br>'
		
		# print (entries_not_found)

		# def generate():
		# 	num = 0
		
		output = makeOutputJson(uniqID, results, mycursor)
		# return app.response_class(generate(), mimetype='text/plain')
		# Save the ignored entries to a json file to be used by datatables
		with open(BASE_DIR+'/webApp/static/predictor/output/'+uniqID+'/ignored.json', 'w') as f:
			json.dump(ignored, f)
		if len(results['predictions']) == 0:
			text = 'None of the variants could not be processed.<br>Please check your input and try again.'
			return render_template('ignored2.html',
								uniqID=json.dumps(uniqID),
								text=json.dumps(text),
								# output=json.dumps(output),
								# error=json.dumps(len(results['entries_not_found']))
								)

		# Save the output to a json file to be used by datatables
		with open(BASE_DIR+'/webApp/static/predictor/output/'+uniqID+'/output.json', 'w') as f:
			json.dump(output, f)
		# Save the results to a json file to be used by the result page
		with open(BASE_DIR+'/webApp/static/predictor/output/'+uniqID+'/results.json', 'w') as f:
			json.dump(results, f)
		return render_template('output2.html',
								uniqID=json.dumps(uniqID),
								output=json.dumps(output),
								error=json.dumps(len(results['entries_not_found']))
								)
		
	@app.route('/result', methods=['GET', 'POST'])
	def result():
		'''
		This function is called by the summary output
		or directly by the user from the browser.
		Takes kinase and mutation as input.
		'''
		if request.method == 'POST':
			data = request.args.get('kinase')
			#print (request.get_json(force=True))
			#return redirect(url_for('home'))
			#return render_template('home.html')
		else:
			
			uniqID = request.args.get('uniqID')
			kinase = request.args.get('kinase')
			mutation = request.args.get('mutation')
			#print (request.args.get('data'))
			results = {}
			if uniqID is None:
				mycursor = connection()
				uniqID = makeUniqID()
				print (uniqID, kinase, mutation)
				results = runPrediction(uniqID, kinase+'/'+mutation)
				output = makeOutputJson(uniqID, results, mycursor)
				# Save the output to a json file to be used by datatables
				with open(BASE_DIR+'/webApp/static/predictor/output/'+uniqID+'/output.json', 'w') as f:
					json.dump(output, f)
			else:
				with open('static/predictor/output/'+uniqID+'/results.json', 'r') as f:
					results = json.load(f)
			if len(results) == 0:
				# return 'No results found' if input is not present in the database
				return render_template('error1.html',
										flaggedInput=json.dumps(kinase+'/'+mutation)
										)
			else:
				return render_template('result2.html',
										uniqID=json.dumps(uniqID),
										kinase=json.dumps(kinase),
										mutation=json.dumps(mutation),
										results=results
										)
	
	@app.route('/ignored/<string:uniqID>', methods=['GET', 'POST'])
	def ignored(uniqID: str):
		'''
		This function is called by the summary output
		or directly by the user from the browser.
		Takes kinase and mutation as input.
		'''
		if request.method == 'POST':
			data = request.args.get('kinase')
			#print (request.get_json(force=True))
			#return redirect(url_for('home'))
			#return render_template('home.html')
		else:
			# uniqID = request.args.get('uniqID')
			# kinase = request.args.get('kinase')
			# mutation = request.args.get('mutation')
			#print (request.args.get('data'))
			# results = {}
			# if uniqID is None:
			# 	uniqID = makeUniqID()
			# 	print (uniqID, kinase, mutation)
			# 	results = runPrediction(uniqID, kinase+'/'+mutation)
			# else:
			# 	with open('static/predictor/output/'+uniqID+'/results.json', 'r') as f:
			# 		results = json.load(f)
			# if len(results) == 0:
			# 	# return 'No results found' if input is not present in the database
			# 	return render_template('error1.html',
			# 							flaggedInput=json.dumps(kinase+'/'+mutation)
			# 							)
			# else:
			with open('static/predictor/output/'+uniqID+'/ignored.json', 'r') as f:
				ignored = json.load(f)
			
			if len(ignored['data']) == 1:
				text = 'The following variants could not be processed.<br>Please check your input and try again.'
			else:
				text = 'The following variants could not be processed.<br>Please check your input and try again.'
			return render_template('ignored2.html',
									uniqID=json.dumps(uniqID),
									text=json.dumps(text),
									# kinase=json.dumps(kinase),
									# mutation=json.dumps(mutation),
									# results=results
									)

	@app.route('/error1')
	def error1():
		return render_template('error1.html')

	@app.route('/AJAXChart', methods=['GET', 'POST'])
	def get_Chart(**kwargs):
		'''
		A function to take uniqID, kinase and mutation as input
		and return prediction information as dic
		'''
		if request.method == 'POST':
			data = request.get_json(force=True)
			uniqID = data['uniqID']
			kinase = data['kinase']
			mutation = data['mutation']
			results = data['results']
		else:
			kinase = kwargs['kinase']
			mutation = kwargs['mutation']
		
		dic = {}
		# print (kinase, mutation)
		activating_prob = results['predictions'][kinase+'/'+mutation]['A']
		deactivating_prob = results['predictions'][kinase+'/'+mutation]['D']
		neutral_prob = results['predictions'][kinase+'/'+mutation]['N']
		resistance_prob = results['predictions'][kinase+'/'+mutation]['RvN']
		activating_AIvLD_prob = results['predictions'][kinase+'/'+mutation]['AIvLD']
		if activating_prob == 'NA':
			activating_AIvLD_prob = 0.0
			deactivating_AIvLD_prob = 0.0
		else:
			deactivating_AIvLD_prob = round(1.0 - float(activating_AIvLD_prob), 3)
		# else:
		# 	activating_AIvLD_prob = float(activating_AIvLD_prob)
		# 	deactivating_AIvLD_prob = 1.0 - activating_AIvLD_prob
		# 	activating_AIvLD_prob = round(activating_AIvLD_prob, 3)
		# 	deactivating_AIvLD_prob = round(deactivating_AIvLD_prob, 3)

		results['predictions'][kinase+'/'+mutation]['activating'] = activating_prob
		results['predictions'][kinase+'/'+mutation]['deactivating'] = deactivating_prob
		results['predictions'][kinase+'/'+mutation]['neutral'] = round(float(neutral_prob), 3)
		results['predictions'][kinase+'/'+mutation]['resistance'] = resistance_prob
		dic = {'activating': results['predictions'][kinase+'/'+mutation]['activating'],
				'deactivating': results['predictions'][kinase+'/'+mutation]['deactivating'],
				'neutral': results['predictions'][kinase+'/'+mutation]['neutral'],
				'resistance': results['predictions'][kinase+'/'+mutation]['resistance'],
				'activating_AIvLD': activating_AIvLD_prob,
				'deactivating_AIvLD': deactivating_AIvLD_prob
				}
		return jsonify(dic)
	
	@app.route('/AJAXSummary', methods=['GET', 'POST'])
	def get_Summary(**kwargs):
		'''
		A function to take uniqID, kinase and mutation as input
		and return summary as dic
		'''
		if request.method == 'POST':
			data = request.get_json(force=True)
			uniqID = data['uniqID']
			kinase = data['kinase']
			mutation = data['mutation']
			results = data['results']
		else:
			kinase = kwargs['kinase']
			mutation = kwargs['mutation']
		with open('static/predictor/output/'+uniqID+'/output.json', 'r') as f:
			output = json.load(f)
		
		text = ''
		for row in output['data']:
			if row['name'] == kinase+'/'+mutation:
				text += '<table><tr><td><b>User input:</b></td><td>'+row['name']+'</td></tr>'
				text += '<tr><tr><td><b>Gene name:</b></td><td>'+row['gene']+'</td></tr>'
				text += '<tr><tr><td><b>UniProt ID:</b></td><td>'+row['uniprot_id']+'</td></tr>'
				# text for acc begins here
				text += '<tr><td><b>UniProt accession:</b></td><td>'
				text += '<a href="https://www.uniprot.org/uniprotkb/'
				text += row['acc']+'/entry" target="_blank">'
				text += row['acc']+'<i class="bi bi-box-arrow-in-up-right"></i></td></tr>'
				# text for acc ends here
				text += '<tr><tr><td><b>Protein name:</b></td><td>'+row['protein_name']+'</td></tr>'
				text += '<tr><td><b>Mutation:</b></td><td>'+row['mutation']+'</td></tr>'
				text += '<tr><td><b>HMM position:</b></td><td>'+row['hmmPos']+'</td></tr>'
				text += '<tr><td><b>Alignment position:</b></td><td>'+row['alnPos']+'</td></tr>'
				text += '<tr><tr><td><b>Region of the site:</b></td><td>'+row['region']+'</td></tr></table>'
				# if row['text'] != '':
				# 	text += '<br><b>More information:</b><br>' + row['text']
				break

		dic = {'text': text}
		return jsonify(dic)
	
	@app.route('/AJAXButtons', methods=['POST', 'GET'])
	def get_Buttons(**kwargs):
		'''
		A function to take uniqID, kinase and mutation as input
		and return summary as dic
		'''
		if request.method == 'POST':
			data = request.get_json(force=True)
			# uniqID = data['uniqID']
			# kinase = data['kinase']
			# mutation = data['mutation']
		else:
			kinase = kwargs['kinase']
			mutation = kwargs['mutation']
		
		buttonText = []
		for line in open('../Create_SVG/Enhancements_May2023/May10th//Ranking_Values.csv', 'r'):
			buttonText.append(line.strip().split(','))

		dic = {'data': buttonText}
		if request.endpoint == 'AJAXButtons':
			return jsonify(dic)
		else:
			return dic
	
	@app.route('/AJAXSummaryTable', methods=['GET', 'POST'])
	def get_SummaryTable(**kwargs):
		'''
		A function to take uniqID, kinase and mutation as input
		and return summary as dic
		'''
		if request.method == 'POST':
			data = request.get_json(force=True)
			uniqID = data['uniqID']
			kinase = data['kinase']
			mutation = data['mutation']
		else:
			kinase = kwargs['kinase']
			mutation = kwargs['mutation']
		with open('static/predictor/output/'+uniqID+'/output.json', 'r') as f:
			output = json.load(f)
		
		text = ''
		for row in output['data']:
			if row['name'] == kinase+'/'+mutation:
				acc = row['acc']
				gene = row['gene']
				hmmPos = row['hmmPos']
				print (acc, gene, hmmPos)
				data, _ = makeText(acc, gene, mutation, hmmPos, connection())
				break
				# text += '<table><tr><td><b>User input:</b></td><td>'+row['name']+'</td></tr>'
				# text += '<tr><tr><td><b>Gene name:</b></td><td>'+row['gene']+'</td></tr>'
				# text += '<tr><tr><td><b>UniProt ID:</b></td><td>'+row['uniprot_id']+'</td></tr>'
				# # text for acc begins here
				# text += '<tr><td><b>UniProt acc:</b></td><td>'
				# text += '<a href="https://www.uniprot.org/uniprotkb/'
				# text += row['acc']+'/entry" target="_blank">'
				# text += row['acc']+'<i class="bi bi-box-arrow-in-up-right"></i></td></tr>'
				# # text for acc ends here
				# text += '<tr><tr><td><b>Protein name:</b></td><td>'+row['protein_name']+'</td></tr>'
				# text += '<tr><td><b>Mutation:</b></td><td>'+row['mutation']+'</td></tr>'
				# text += '<tr><tr><td><b>Region of the site:</b></td><td>'+row['region']+'</td></tr></table>'
				# if row['text'] != '':
				# 	text += '<br><b>More information:</b><br>' + row['text']
				# break

		dic = {'data': data, 'acc': acc, 'gene': gene}
		return jsonify(dic)
	
	@app.route('/AJAXAlignment', methods=['GET', 'POST'])
	def get_Alignment(**kwargs):
		'''
		A function to take uniqID, kinase and mutation as input
		and return summary as dic
		'''
		if request.method == 'POST':
			data = request.get_json(force=True)
			uniqID = data['uniqID']
			kinase = data['kinase']
			mutation = data['mutation']
			results = data['results']
			ws = data['WS']
			topN = data['topN']
			sortTypeText = data['sortTypeText']
		else:
			kinase = kwargs['kinase']
			mutation = kwargs['mutation']

		with open('static/predictor/output/'+uniqID+'/output.json', 'r') as f:
			output = json.load(f)
		
		accs_in_alignment = {}
		for line in open('static/hmm/humanKinasesTrimmed.clustal', 'r'):
			if line.startswith('sp'):
				accs_in_alignment[line.split('|')[1]] = line.split()[0]

		mycursor = connection(database='kinase_project2')

		dic_mutations_info = {}
		# mut_type_name = {'A': 'Activating', 'D': 'Deactivating', 'R': 'Resistance'}
		mut_type_name = {'activating': 'Activating', 'increase': 'Activating', \
		   				'loss': 'Deactivating', 'decrease': 'Deactivating',\
						'resistance': 'Resistance'}
		ptm_type_name = {'me': 'Methylation', 'm1': 'Methylation', 'm2': 'Methylation',\
		   				'm3': 'Methylation', 'p': 'Phosphorylation', 'ac': 'Acetylation',
						'ub': 'Ubiquitination', 'sm': 'Sumoylation', 'gl': 'O-GlcNAc'}
		for row in output['data']:
			if row['name'] != kinase+'/'+mutation: continue			
			mutation_position = int(row['mutation'][1:-1])
			## Mutations
			mycursor.execute("SELECT mutation, wtpos, mut_type, acc, gene, info, pubmed FROM mutations")
			hits = mycursor.fetchall()
			for hit in hits:
				otherMutation, position, mutType, acc, gene, info, pubmedIDs = hit
				if acc not in dic_mutations_info:
					# dic_mutations_info[acc] = {'A':[], 'D':[], 'R':[]}
					dic_mutations_info[acc] = {'Activating':[], 'Deactivating':[], 'Resistance':[]}
				if mutType == 'neutral': continue
				if acc not in accs_in_alignment: continue
				if mutType in ['activating', 'increase', 'loss', 'decrease']:
					# pubmed_ids = extract_pubmed_ids(info.replace('"', '')) # remove double quotes
					pubmed_ids = pubmedIDs.split(',')
					text = ''
					if pubmed_ids != []:
						search_text = '+or+'.join(pubmed_ids)
						# text += '<a href=\"https://pubmed.ncbi.nlm.nih.gov/?term=' + str(search_text) + '\" target=\"_blank\">PubMed<i class="bi bi-box-arrow-in-up-right"></i></a>'
						text += 'https://pubmed.ncbi.nlm.nih.gov/?term=' + str(search_text)
				else:
					# text = ' ' + '<a href=\"https://cancer.sanger.ac.uk/cosmic/gene/analysis?ln='
					# text += gene +'#drug-resistance\" target=\"_blank\">COSMIC <i class="bi bi-box-arrow-in-up-right"></i></a>'
					text = 'https://cancer.sanger.ac.uk/cosmic/gene/analysis?ln='+gene +'#drug-resistance'
				name = accs_in_alignment[acc]
				start = int(name.split('|')[-1])
				if int(position) >= start:
					# dic_mutations_info[acc][mut_type_name[mutType]].append(str(position))
					dic_mutations_info[acc][mut_type_name[mutType]].append([str(position), otherMutation, text])
			# PTMs
			mycursor.execute("SELECT uniprotaa, uniprotpos, ptmtype, acc FROM ptms")
			hits = mycursor.fetchall()
			for hit in hits:
				wtAA, position, ptmType, acc = hit
				if acc not in accs_in_alignment: continue
				if acc not in dic_mutations_info:
					dic_mutations_info[acc] = {}
				name = accs_in_alignment[acc]
				start = int(name.split('|')[-1])
				if int(position) >= start:
					if ptm_type_name[ptmType] not in dic_mutations_info[acc]: dic_mutations_info[acc][ptm_type_name[ptmType]] = []
					# text = '<a href=\"http://www.phosphosite.org/uniprotAccAction?id='+ acc +'\" target=\"_blank\">PhosphoSitePlus <i class="bi bi-box-arrow-in-up-right"></i></a>'
					text = 'http://www.phosphosite.org/uniprotAccAction?id='+ acc
					dic_mutations_info[acc][ptm_type_name[ptmType]].append([str(position), wtAA+str(position), text])
					# dic_mutations_info[acc][ptm_type_name[ptmType]].append(str(position))
			break
		
		alignment = 'static/hmm/humanKinasesTrimmed.clustal'
		# alignment = 'static/hmm/humanKinasesHitsSplitTrimmedWeb.aln'
		# alignment = BASE_DIR + '/alignments/humanKinasesHitsSplitTrimmedWeb.aln'
		alignment = BASE_DIR + '/alignments/humanKinasesPkinasePK_Tyr_Ser-ThrAll_no_gapsTrimmedWeb.aln'
		# resetDic(dic_mutations_info, alignment)
		dic_mutations_info = resetDic(dic_mutations_info, alignment)
		# print (dic_mutations_info)
		# print (row['acc'], mutation_position)
		'''
		filename = create_svg_20230428_kinases_GS.main('static/hmm/humanKinasesTrimmed.clustal',\
					row['acc'], mutation_position, int(ws), int(topN), dic_mutations_info, \
						path = 'static/predictor/output/'+uniqID+'/')
		'''
		with open(conservation_dic_path) as g:
			overconserv = g.read()
		overallconservation = ast.literal_eval(overconserv)
		with open(identity_dic_path) as h:
			data_ident = h.read()
		identitydictionary = ast.literal_eval(data_ident)
		sortingvalue = '1'
		dic_buttons = get_Buttons()
		for values in dic_buttons['data']:
			if values[1] == sortTypeText:
				sortingvalue = str(values[0])
				break
		print (f'sortingValue is {sortingvalue}')
		# geeky_file = open('sample_dic_mutation_info.txt', 'wt')
		# geeky_file.write(str(dic_mutations_info))
		# try:

		# The following searching  is needed to find the name of entry
		# in the alignment file that contains the mutation position
		# entry is like GN|Acc|start-end
		entry_to_search = ''
		for entry in dic_mutations_info:
			entry_acc = entry.split('|')[1]
			start, end = entry.split('|')[-1].split('-')
			if int(start) <= int(mutation_position) <= int(end) and entry_acc == row['acc']:
				entry_to_search = entry
				break
		feature_dic = {}
		for line in gzip.open('../alignments/humanKinasesPkinasePK_Tyr_Ser-ThrAll_no_gapsTrimmed_ss.tsv.gz', 'rt', encoding='utf-8'):
			if line.startswith('#'): continue
			name = str(line.split()[0])
			start, end = line.split()[1].rstrip().split('-')
			# print (name, start, end)
			feature_dic[name] = [i for i in range(int(start), int(end)+1)]
		try:
			filename = create_svg.main(sortingvalue, identitydictionary,\
					overallconservation, alignment, entry_to_search,\
					mutation_position, int(ws), int(topN), dic_mutations_info,\
					feature_dic,
					path = 'static/predictor/output/'+uniqID+'/')
		except Exception as e:
			print (e)
			filename = ''

		runStatus = 'success' if filename != '' else 'error'
		print (runStatus)

		dic = {
			'filepath': 'static/predictor/output/'+uniqID+'/'+filename,
	 		'status': runStatus
			}
		return jsonify(dic)

	@app.route('/AJAXModals', methods=['GET', 'POST'])
	def AJAXModals (**kwargs):
		if request.method == 'POST':
			data = request.get_json(force=True)
			modalTitle = data['modalTitle']
		else:
			kinase = kwargs['kinase']
			mutation = kwargs['mutation']
		dic = {}
		modalText = ''
		for line in open('static/modals/'+modalTitle+'Modal.txt', 'r'):
			modalText += line
		dic['modalText'] = modalText
		print (modalText)
		return jsonify(dic)
	
	@app.route('/AJAXHomology', methods=['GET', 'POST'])
	def get_homology (**kwargs):
		'''
		mycursor = connection()
		mycursor.execute("SELECT mutation, mut_type, acc FROM mutations")
		hits = mycursor.fetchall()
		dic_mutations_info = {'A':[], 'D':[], 'R':[], 'N':[]}
		for hit in hits:
			mutation, mutType, acc = hit
			dic_mutations_info[mutType].append(acc+'/'+mutation)
		
		# extract homology scores from DB
		homologs = ['orth', 'all_homs', 'bpsh', 'bpso', 'excl_para', 'spec_para']
		# dic_scores = {'A':[], 'D':[], 'R':[], 'N':[]}
		dic_scores = {}
		data = []
		for homolog in homologs:
			if homolog not in dic_scores: dic_scores[homolog] = {'A':[], 'D':[], 'R':[], 'N':[]}
			for mutType in dic_mutations_info:
				# if mutType not in ['A', 'D']: continue
				for instance in dic_mutations_info[mutType]:
					acc = instance.split('/')[0]
					mutation = instance.split('/')[1]
					position = str(mutation[1:-1])
					wtAA = mutation[0]
					mutAA = mutation[-1]
					mycursor.execute("SELECT "+wtAA+"_score, "+mutAA+"_score  FROM "+homolog+"\
									WHERE acc = '"+acc+"' AND position = '"+position+"'")
					hits = mycursor.fetchall()
					for hit in hits:
						wt_score, mut_score = hit
						dic_scores[homolog][mutType].append(float(mut_score) - float(wt_score))
						data.append([homolog, mutType, acc, mutation, str(float(mut_score)-float(wt_score))])
		
		text = 'homolog\tmutType\tacc\tmutation\tscore\n'
		for row in data:
			text += '\t'.join(row) + '\n'
		open('static/data/homologyScores.tsv', 'w').write(text)
					
			
		dic = dic_scores
		'''
		dic = {}
		return jsonify(dic)
	
	@app.route('/AJAXPTM', methods=['GET', 'POST'])
	def get_PTM (**kwargs):
		'''
		PTM_TYPES = ['ac', 'gl', 'm1', 'm2', 'm3', 'me', 'p', 'sm', 'ub']
		mycursor = connection()
		mycursor.execute("SELECT mutation, mut_type, acc FROM mutations")
		hits = mycursor.fetchall()
		dic_mutations_info = {'A':[], 'D':[], 'R':[], 'N':[]}
		for hit in hits:
			mutation, mutType, acc = hit
			dic_mutations_info[mutType].append(acc+'/'+mutation)
		
		WS = 3
		if WS > 0: ws = WS - 1
		ws = int(WS/2)
		# extract homology scores from DB
		ptm_header = []
		for position in range(ws*(-1), ws+1):
			for ptm_type in PTM_TYPES:
				ptm_header.append(ptm_type+'_'+str(position))
				ptm_header.append(ptm_type+'_pfam_'+str(position))
		
		dic_scores = {}
		data = []
		for mutType in dic_mutations_info:
			# if mutType not in ['A', 'D']: continue
			for instance in dic_mutations_info[mutType]:
				acc = instance.split('/')[0]
				mutation = instance.split('/')[1]
				position = int(mutation[1:-1])
				wtAA = mutation[0]
				mutAA = mutation[-1]
				values = fetchData.getPTMscore(mycursor, acc, position, WS)
				for value, ptm_head in zip(values, ptm_header):
					if ptm_head not in dic_scores: dic_scores[ptm_head] = {'A':[], 'D':[], 'R':[], 'N':[]}
					dic_scores[ptm_head][mutType].append(value)
					data.append([ptm_head, mutType, acc, mutation, str(value)])
		
		text = 'ptm_head\tmutType\tacc\tmutation\tscore\n'
		for row in data:
			text += '\t'.join(row) + '\n'
		open('static/data/ptmScores.tsv', 'w').write(text)
		'''
			
		# dic = dic_scores
		
		dic = {}
		return jsonify(dic)
	
	@app.route('/AJAXADR', methods=['GET', 'POST'])
	def get_ADR (**kwargs):
		'''
		MUT_TYPES = ['A', 'D', 'R']
		mycursor = connection()
		mycursor.execute("SELECT mutation, mut_type, acc FROM mutations")
		hits = mycursor.fetchall()
		dic_mutations_info = {'A':[], 'D':[], 'R':[], 'N':[]}
		for hit in hits:
			mutation, mutType, acc = hit
			dic_mutations_info[mutType].append(acc+'/'+mutation)
		
		WS = 3
		if WS > 0: ws = WS - 1
		ws = int(WS/2)
		# extract homology scores from DB
		mut_header = []
		for position in range(ws*(-1), ws+1):
			for mut_type in MUT_TYPES:
				mut_header.append(mut_type+'_'+str(position))
				mut_header.append(mut_type+'_pfam_'+str(position))
		
		dic_scores = {}
		data = []
		num = 0
		for mutType in dic_mutations_info:
			# print (num)
			# if mutType not in ['A', 'D']: continue
			for instance in dic_mutations_info[mutType]:
				acc = instance.split('/')[0]
				mutation = instance.split('/')[1]
				position = int(mutation[1:-1])
				wtAA = mutation[0]
				mutAA = mutation[-1]
				values = fetchData.getADRvector(mycursor, acc, position, WS, WS)
				print (values)
				for value, mut_head in zip(values, mut_header):
					if mut_head not in dic_scores: dic_scores[mut_head] = {'A':[], 'D':[], 'R':[], 'N':[]}
					dic_scores[mut_head][mutType].append(value)
					data.append([mut_head, mutType, acc, mutation, str(value)])
		
		text = 'adr_head\tmutType\tacc\tmutation\tscore\n'
		for row in data:
			text += '\t'.join(row) + '\n'
		open('static/data/adrScores.tsv', 'w').write(text)
		'''
			
		# dic = dic_scores
		
		dic = {}
		return jsonify(dic)

	@app.route('/progress')
	def progress():
		return render_template('progress.html')
	
	@app.route('/progress2', methods=['GET', 'POST'])
	def progress2():
		uniqID = makeUniqID()
		return render_template('progress2.html', uniqID=uniqID)
	
	@app.route('/test')
	def test():
		return render_template('test.html')

'''
This will be called if you run this from command line
'''
if __name__ == "__main__":
	app = Flask(__name__)
	configureRoutes(app)
	app.run(debug=True)
	print ('Ciao!')
else:
	#import webApp.callUniProtAPI as callUniProtAPI
	app = Flask(__name__)
	configureRoutes(app)
