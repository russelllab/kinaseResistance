#!/usr/bin/env python3.10

'''
A script to test the functions in prepareTestData.py
'''

import pytest
import os, sys
os.sys.path.append('../ML/')
import prepareTestData

def test_predict():
    """
    Test the predict function in prepareTestData.py
    """
    results = prepareTestData.predict('sample_mutations3.txt', BASE_DIR = '../')
    assert 'BRAF/T600E' in results['entries_not_found']
    assert '-' == results['predictions']['BRAF/V600E']['ptmType']
    assert 'p' == results['predictions']['BRAF/T599E']['ptmType']
    assert 'NA' == results['predictions']['FGFR2/C278F']['prediction']
    assert float(results['predictions']['P46734/S218E']['prediction']) > 0.5

if __name__ == '__main__':
    test_predict()