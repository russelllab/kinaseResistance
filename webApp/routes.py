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
	BASE_URL = 'http://kinaser.russelllab.org/'
	# BASE_DIR = '/net/home.isilon/ds-russell/kinaseResistance/'
	BASE_DIR = '/var/www/flask_apps/kinaseResistance/'
else:
	BASE_URL = 'http://127.0.0.1:5000/'
	BASE_DIR = '../'
	BASE_DIR = '/home/gurdeep/projects/kinaseResistance/'


sys.path.insert(1, BASE_DIR+'/ML/')
import prepareTestData

def connection():
    '''Function to connect to postgresql database'''
    mydb = psycopg2.connect(
                            database = "kinase_project",
                            user = "gurdeep",
                            password = "hellokitty",
                            host = "localhost",
                            port = "5432")
    return mydb

def makeUniqID():
	'''
	Generates a unique ID for the job
	'''
	# A unique ID is created before the pocess begins
	# Each submitted job is assigned a unique ID
	uniqID = ''.join(random.choices(string.ascii_uppercase + string.digits, k=5))
	return uniqID

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

	@app.route('/test')
	def test():
		'''
		Test
		'''
		return 'Hello World'

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
		return render_template('home.html', uniqID=uniqID)

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

		mydb = connection()
		mydb.autocommit = True
		mycursor = mydb.cursor()
		#sys.exit()
		# inputMuts = request.form['inputMut']
		results = runPrediction(uniqID, inputMuts)
		# print (results)
		if isinstance(results, float):
			return app.response_class(results, mimetype='text/plain')

		if len(results) == 0:
			return render_template('error1.html',
								flaggedInput=json.dumps(kinase+'/'+mutation)
								)
		
		output = []
		entries_not_found = results['entries_not_found']
		error_html_text = ''
		for name in entries_not_found:
			error_html_text += name + ' ' + entries_not_found[name] + '<br>'
		
		# print (entries_not_found)

		# def generate():
		# 	num = 0
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
			dic['view'] += '\" target=\"_blank\">View</a>'
			dic['gene'] = results['predictions'][name]['gene']
			dic['ptmType'] = results['predictions'][name]['ptmType']
			dic['mutType'] = results['predictions'][name]['mutType']
			dic['prediction'] = results['predictions'][name]['prediction']
			dic['hmmPos'] = results['predictions'][name]['hmmPos']
			output.append(dic)
			# yield str(num) + '\n'
		
		# return app.response_class(generate(), mimetype='text/plain')
	
		output = {'data': output}
		# Save the output to a json file to be used by datatables
		with open(BASE_DIR+'/webApp/static/predictor/output/'+uniqID+'/output.json', 'w') as f:
			json.dump(output, f)
		# Save the results to a json file to be used by the result page
		with open(BASE_DIR+'/webApp/static/predictor/output/'+uniqID+'/results.json', 'w') as f:
			json.dump(results, f)
		return render_template('output.html',
								uniqID=json.dumps(uniqID),
								output=json.dumps(output),
								error=json.dumps(error_html_text)
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
				uniqID = makeUniqID()
				print (uniqID, kinase, mutation)
				results = runPrediction(uniqID, kinase+'/'+mutation)
			else:
				with open('static/predictor/output/'+uniqID+'/results.json', 'r') as f:
					results = json.load(f)
			if len(results) == 0:
				# return 'No results found' if input is not present in the database
				return render_template('error1.html',
										flaggedInput=json.dumps(kinase+'/'+mutation)
										)
			else:
				return render_template('result.html',
										uniqID=json.dumps(uniqID),
										kinase=json.dumps(kinase),
										mutation=json.dumps(mutation),
										results=results
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
		activating_prob = results['predictions'][kinase+'/'+mutation]['prediction']
		if activating_prob == 'NA':
			activating_prob = 0.0
			deactivating_prob = 0.0
		else:
			activating_prob = float(activating_prob)
			deactivating_prob = 1.0 - activating_prob

		results['predictions'][kinase+'/'+mutation]['activating'] = activating_prob
		results['predictions'][kinase+'/'+mutation]['deactivating'] = deactivating_prob
		results['predictions'][kinase+'/'+mutation]['neutral'] = 0.5
		results['predictions'][kinase+'/'+mutation]['resistant'] = 0.5
		dic = {'activating': results['predictions'][kinase+'/'+mutation]['activating'],
				'deactivating': results['predictions'][kinase+'/'+mutation]['deactivating'],
				'neutral': results['predictions'][kinase+'/'+mutation]['neutral'],
				'resistant': results['predictions'][kinase+'/'+mutation]['resistant']}
		return jsonify(dic)
		
	@app.route('/progress')
	def progress():
		return render_template('progress.html')

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
