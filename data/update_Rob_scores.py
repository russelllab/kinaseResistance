import os, gzip

path = '../KA/robScores2/'
for files in os.listdir(path):
    l = ''
    for line in gzip.open(path+files, 'rt'):
        if line[0] == '#':
            l += line.replace('\n','') + '\tUniProt_Acc' + '\n'
        else:
            l += line.replace('\n','') + '\t' + files.split('.')[0] + '\n'
    gzip.open('../KA/robScores3/'+files, 'wt').write(l)