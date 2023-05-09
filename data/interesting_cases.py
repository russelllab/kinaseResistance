import os
import sys
sys.path.insert(0, ('../ML/'))
import fetchData

mydb = fetchData.connection()
mydb.autocommit = True
mycursor = mydb.cursor()

mycursor.execute("SELECT gene, mut_type, mutation, wtpos, info FROM mutations where wtaa = 'S' or wtaa = 'T'")
hits = mycursor.fetchall()
dict_hits = {}
dic_info = {}
for hit in hits:
    gene = hit[0]
    mut_type = hit[1]
    if mut_type not in ['A', 'D']: continue
    mutation = hit[2]
    wtpos = hit[3]
    info = hit[4]
    # dic[gene + '/' + mutation] = info
    name = gene + '/' + mutation[:-1]
    if name not in dict_hits:
        dict_hits[name] = {'A': [], 'D': []}
    dict_hits[name][mut_type].append(mutation)

for name in dict_hits:
    if dict_hits[name]['A'] != [] and dict_hits[name]['D'] != []:
        print(name, dict_hits[name])