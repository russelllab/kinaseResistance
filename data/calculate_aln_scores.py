'''
calculate log-odds scores using Rob's method
'''
import os, gzip
acc = 'P15056'
aln_category = 'orth'
path_to_aln = '../alignments/' + acc[:4] + '/'
pfam_dom = 'Pkinase'
pfam_path = '../pfam/'

'''
Extract FASTA from given MSA
'''
aln_dic = {}
for files in os.listdir(path_to_aln):
    if files.split('_')[-1] != aln_category+'.aln.gz':
        continue
    if files.split('_')[0] != acc:
        continue
    for line in gzip.open(path_to_aln+files, 'rt'):
        if len(line.split()) <= 1:
            continue
        elif line.split()[0] == 'CLUSTAL':
            continue
        name = line.split()[0]
        msa_seq = line.split()[1].rstrip('\n')
        if name not in aln_dic: aln_dic[name] = ''
        aln_dic[name] += msa_seq.replace('-', '')

l = ''
for name in aln_dic:
    # print (name, aln_dic[name])
    l += '>' + name + '\n'
    l += aln_dic[name] + '\n'

gzip.open(path_to_aln + acc + '_' + aln_category + '.fasta.gz', 'wt').write(l)

'''
Perform HMMalign
'''

os.system('zcat ' + path_to_aln + acc + '_' + aln_category + '.fasta.gz|'
        'hmmalign -o '+
        path_to_aln + acc + '_' + aln_category + '.hmmaln ' + 
        pfam_path+pfam_dom+'.hmm -')