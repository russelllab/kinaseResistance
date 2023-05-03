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
import re, sys, math
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
	BASE_URL = 'http://kinasex.russelllab.org/'
	# BASE_DIR = '/net/home.isilon/ds-russell/kinaseResistance/'
	BASE_DIR = '/var/www/flask_apps/kinaseResistance/'
else:
	BASE_URL = 'http://127.0.0.1:5000/'
	BASE_DIR = '../'
	BASE_DIR = '/home/gurdeep/projects/kinaseResistance/'


sys.path.insert(1, BASE_DIR+'/ML/')
import prepareTestData

sys.path.insert(1, BASE_DIR+'/Create_SVG/Vlatest/')
# import create_svg_20230426_kinases_GS
# import create_svg_20230428_kinases_GS
import create_svg_20230502_kinases_GS

def connection():
    '''Function to connect to postgresql database'''
    mydb = psycopg2.connect(
                            database = "kinase_project",
                            user = "gurdeep",
                            password = "hellokitty",
                            host = "localhost",
                            port = "5432")
    
    mydb.autocommit = True
    mycursor = mydb.cursor()
    return mycursor

def makeText(acc, mutation, mycursor):
	'''
	Make text for prediction
	'''
	ws = 3
	if ws > 0: ws -= 1
	ws = int(ws/2)
	dic_mutations = {'A': 'activating', 'D': 'deactivating', 'R': 'resistance'}
	dic_ptms = {'p': 'phosphorylation', 'ub': 'ubiquitination', 'ac': 'acetylation', 'me': 'methylation', 'gl': 'glycosylation', 'sm': 'sumoylation', 'm1': 'myristoylation', 'm2': 'palmitoylation', 'm3': 'myristoylation'}
	text = ''
	mutation_position = int(mutation[1:-1])
	for position in range(mutation_position-ws, mutation_position+ws+1):
		mycursor.execute("SELECT uniprotaa, ptmtype FROM ptms \
						WHERE acc = %s and uniprotpos = %s", (acc, str(position)))
		hits = mycursor.fetchall()
		for entry in hits:
			text += "<b>"+str(entry[0]) + str(position)+ "</b>" + ' is a known ' + dic_ptms[entry[1]] + ' site<br>'
	
	for position in range(mutation_position-ws, mutation_position+ws+1):
		mycursor.execute("SELECT mutation, mut_type, info FROM mutations \
						WHERE acc = %s and wtpos = %s", (acc, str(position)))
		hits = mycursor.fetchall()
		for entry in hits:
			text += "<b>" + str(entry[0]) + "</b>" + ' is a known '+dic_mutations[entry[1]]+' mutation.'
			text += ' Ref: ' + entry[2] + '<br>'
	
	return text

def makeUniqID():
	'''
	Generates a unique ID for the job
	'''
	# A unique ID is created before the pocess begins
	# Each submitted job is assigned a unique ID
	uniqID = ''.join(random.choices(string.ascii_uppercase + string.digits, k=5))
	return uniqID

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
		dic['acc'] = results['predictions'][name]['acc']
		dic['view'] = '<a href=\"'
		dic['view'] += BASE_URL
		dic['view'] += 'result?uniqID='+uniqID+'&kinase=' + kinase + '&mutation=' + mutation
		# dic['view'] += '\" target=\"_blank\">View</a>'
		dic['view'] += '\" target=\"_blank\"><i class="bi bi-box-arrow-in-up-right"></i></a>'
		dic['gene'] = results['predictions'][name]['gene']
		dic['ptmType'] = results['predictions'][name]['ptmType']
		dic['mutType'] = results['predictions'][name]['mutType']
		dic['predAD'] = results['predictions'][name]['predAD']
		dic['predRN'] = results['predictions'][name]['predRN']
		dic['hmmPos'] = results['predictions'][name]['hmmPos']
		dic['text'] = makeText(dic['acc'], dic['mutation'], mycursor)
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
	results = prepareTestData.predict(BASE_DIR+'/webApp/static/predictor/output/'+uniqID+'/input.txt', \
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
		uniqID = makeUniqID()
		# return 'Home'
		# return render_template('home.html', uniqID=uniqID)
		return render_template('progress2.html', uniqID=uniqID)

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
			f = open(BASE_DIR+'/tests/sample_mutations3.txt', "r")
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
			kinase, mutation = name.split('/')
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
			text = 'Unfortunately, none of the mutations you submitted were found in the database.<br>Please check the input and try again.'
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
				text = 'Unfortunately, the following mutation you submitted was not found in the database.<br>Please check the input and try again.'
			else:
				text = 'Unfortunately, the following mutations you submitted were not found in the database.<br>Please check the input and try again.'
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
		activating_prob = results['predictions'][kinase+'/'+mutation]['predAD']
		if activating_prob == 'NA':
			activating_prob = 0.0
			deactivating_prob = 0.0
		else:
			activating_prob = float(activating_prob)
			deactivating_prob = 1.0 - activating_prob

		resistance_prob = results['predictions'][kinase+'/'+mutation]['predRN']
		results['predictions'][kinase+'/'+mutation]['activating'] = activating_prob
		results['predictions'][kinase+'/'+mutation]['deactivating'] = deactivating_prob
		results['predictions'][kinase+'/'+mutation]['neutral'] = 0.5
		results['predictions'][kinase+'/'+mutation]['resistant'] = resistance_prob
		dic = {'activating': results['predictions'][kinase+'/'+mutation]['activating'],
				'deactivating': results['predictions'][kinase+'/'+mutation]['deactivating'],
				'neutral': results['predictions'][kinase+'/'+mutation]['neutral'],
				'resistant': results['predictions'][kinase+'/'+mutation]['resistant']}
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
				text += '<table><tr><td><b>Input:</b></td><td>'+row['name']+'</td></tr>'
				text += '<tr><tr><td><b>Gene:</b></td><td>'+row['gene']+'</td></tr>'
				text += '<tr><td><b>UniProt Acc:</b></td><td>'+'<a href="https://www.uniprot.org/uniprotkb/'+row['acc']+'/entry" target="_blank">'+row['acc']+'</td></tr>'
				text += '<tr><td><b>Mutations:</b></td><td>'+row['mutation']+'</td></tr></table>'
				if row['text'] != '':
					text += '<br><b>More information:</b><br>' + row['text']
				break

		dic = {'text': text}
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
		else:
			kinase = kwargs['kinase']
			mutation = kwargs['mutation']

		with open('static/predictor/output/'+uniqID+'/output.json', 'r') as f:
			output = json.load(f)
		
		accs_in_alignment = {}
		for line in open('static/hmm/humanKinasesTrimmed.clustal', 'r'):
			if line.startswith('sp'):
				accs_in_alignment[line.split('|')[1]] = line.split()[0]

		mycursor = connection()

		text = ''
		for row in output['data']:
			if row['name'] != kinase+'/'+mutation: continue			
			mutation_position = int(row['mutation'][1:-1])
			mycursor.execute("SELECT wtpos, mut_type, acc FROM mutations")
			hits = mycursor.fetchall()
			dic_mutations_info = {}
			for hit in hits:
				position, mutType, acc = hit
				if acc not in dic_mutations_info:
					dic_mutations_info[acc] = {'A':[], 'D':[], 'R':[]}
				if mutType == 'N': continue
				if acc not in accs_in_alignment: continue
				name = accs_in_alignment[acc]
				start = int(name.split('|')[-1])
				# start, end = name.split('|')[-1].split('-')
				# if start == 'start': start = 1
				# else: start = int(start)
				# if end == 'end': end = 10000
				# else: end = int(end)
				# if int(position) >= start and int(position) <= end:
				if int(position) >= start:
					dic_mutations_info[acc][mutType].append(str(position))
			break
		
		# print (dic_mutations_info)
		# print (row['acc'], mutation_position)
		'''
		filename = create_svg_20230428_kinases_GS.main('static/hmm/humanKinasesTrimmed.clustal',\
					row['acc'], mutation_position, int(ws), int(topN), dic_mutations_info, \
						path = 'static/predictor/output/'+uniqID+'/')
		'''
		try:
			filename = create_svg_20230502_kinases_GS.main('static/hmm/humanKinasesTrimmed.clustal',\
					row['acc'], mutation_position, int(ws), int(topN), dic_mutations_info, \
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

	@app.route('/progress')
	def progress():
		return render_template('progress.html')
	
	@app.route('/progress2')
	def progress2():
		return render_template('progress2.html')
	
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
