#!/usr/bin/env python3.10

'''
A script to test the functions in prepareTestData.py
'''

import pytest
import os, sys
import pandas as pd
os.sys.path.append('ML/')
import prepareTestData
os.sys.path.append('webApp/')
from routes import configureRoutes, makeUniqID
from routes import *
from flask import Flask

@pytest.fixture
def client():
    """
    Define a fixture for the Flask test client
    """
    app = Flask(__name__, template_folder='../webApp/templates/')
    ## Call routes
    configureRoutes(app)
    ## Create a client
    client = app.test_client()

    return client

def test_index(client):
    """
    Test that the index page returns a 200 status code
    """
    client = client
    response = client.get('/')
    print (response)
    assert response.status_code == 200
    assert b'KinAct' in response.data

def test_output(client):
    """
    Test that the output page returns a 200 status code
    and has BRAF/V600E in the output
    """
    uniqID = makeUniqID()
    print ('---'+uniqID+'---')
    client = client
    response = client.get('/output/'+uniqID)
    # print (response.data)
    # print (response.data.decode('utf-8').split('<table')[1])
    assert response.status_code == 200
    assert b'BRAF/V600E' in response.data
