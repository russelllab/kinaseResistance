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
	#BASE_DIR = '/var/www/flask_apps/mechismoX/'
	BASE_DIR = '/net/home.isilon/ds-russell/kinaseResistance/'
else:
	BASE_URL = 'http://127.0.0.1:5000/'
	BASE_DIR = '../'


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
    
def retrieve_entries(protein, mutation, mycursor):
	'''For a give acc/gene/uniprot ID, retireve known information'''
	mycursor.execute("select acc, gene from kinases where acc=%s", (protein,))
	hits = mycursor.fetchone()
	if hits is None:
		mycursor.execute("select acc, gene from kinases where gene=%s", (protein,))
		hits = mycursor.fetchone()
		if hits is None:
			print (f'Neither acc nor gene with name {protein} found')
			return None, None, None, None
	acc, gene = hits[0], hits[1]
	
	## fetch ptm_types
	mycursor.execute(\
					'select ptmtype from ptms \
					where acc=%s and uniprotpos=%s', \
					(acc, mutation[1:-1],)\
					)
	hits = mycursor.fetchone()
	if hits is not None: ptmType = hits[0]
	else: ptmType = 'None'
	
	## fetch mut_types
	mycursor.execute(\
					'select mut_type from mutations \
					where acc=%s and mutation=%s', \
					(acc, mutation,)\
					)
	hits = mycursor.fetchone()
	if hits is not None: mutType = hits[0]
	else: mutType = 'None'

	print (acc, gene, ptmType, mutType)
	return (acc, gene, ptmType, mutType)

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
		uniqID = ''.join(random.choices(string.ascii_uppercase + string.digits, k=5))
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
			mydb = connection()
			mydb.autocommit = True
			mycursor = mydb.cursor()
			#sys.exit()
			inputMuts = request.form['inputMut']
			if os.path.isfile('static/predictor/output/'+uniqID) is False:
				os.system('mkdir static/predictor/output/'+uniqID)
			with open('static/predictor/output/'+uniqID+'/input.txt', 'w') as f:
				f.write(inputMuts)
			results = prepareTestData.predict('static/predictor/output/'+uniqID+'/input.txt', \
				     BASE_DIR = BASE_DIR)
			# return inputMuts
			output = []; seq = ''; dic = {}
			for line in inputMuts.split('\n'):
				if line[0] == '#' or line.lstrip().rstrip() == '':
					continue
				protein = line.split('/')[0].lstrip()
				mutation = line.split('/')[1].rstrip()
				acc, gene, ptmType, mutType = retrieve_entries(protein, mutation, mycursor)
				if acc is None: continue
				name = acc + '/' + mutation
				dic = {}
				dic['mutation'] = mutation
				dic['acc'] = acc
				dic['gene'] = gene
				dic['ptmType'] = ptmType
				dic['mutType'] = mutType
				if name in results:
					dic['prediction'] = results[name]['prediction']
					dic['hmmPos'] = results[name]['hmmPos']
				else:
					dic['prediction'] = '-'
					dic['hmmPos'] = '-'
				output.append(dic)
			output = {'data': output}
			with open('static/predictor/output/'+uniqID+'/output.json', 'w') as f:
				json.dump(output, f)
			return render_template('output.html', uniqID=json.dumps(uniqID), output=json.dumps(output))
		else:
			print ('GET')
			return 'You called GET'

'''
This will be called if you run this from command line
'''
if __name__ == "__main__":
	# import callUniProtAPI
	app = Flask(__name__)
	configureRoutes(app)
	app.run(debug=True)
	print ('Ciao!')
else:
	#import webApp.callUniProtAPI as callUniProtAPI
	# import callUniProtAPI
	app = Flask(__name__)
	configureRoutes(app)
	# app.run(debug=True)
