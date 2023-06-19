#!/usr/bin/env python3
# coding: utf-8

import argparse
import os
import re

parser = argparse.ArgumentParser(description='Prepare alignments')
parser.add_argument('-i', '--input', help='input file', required=True)

args = parser.parse_args()
input_file = args.input

dic = {}
for line in open(input_file, 'r'):
    if line.startswith('#') or line.startswith('//'): continue
    if line.lstrip().rstrip() == '': continue
    # print (line.split())
    name = line.split()[0]
    seq = line.split()[1].rstrip().lstrip()
    if name not in dic: dic[name] = ''
    dic[name] += seq

for name in dic:
    print ('>' + name)
    print (dic[name])