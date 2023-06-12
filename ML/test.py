#!/usr/bin/env python3
# coding: utf-8

import fetchData

mydb = fetchData.connection('kinase_project2')
mydb.autocommit = True
mycursor = mydb.cursor()

kinases = {}
acc = 'P11309'
mutation_position = 97
WS = 5
adr_row = fetchData.getADRvector(mycursor, acc, mutation_position, kinases, ws=WS)
adr_row = fetchData.getCountAAchange(mycursor, acc, mutation_position, kinases, ws=WS)
print (adr_row)
AA = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',\
    'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
for i in range(WS*3):
    for j in range(20):
        print (i-2, AA[j], adr_row[20*i+j])
# print (adr_row)