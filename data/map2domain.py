'''
Script to map UniProt seq positions
to Pkinase domain positions
'''
import os, gzip, sys
import argparse, threading
import os, sys, gzip

parser = argparse.ArgumentParser(
                    prog='MapUniProt2HMM',
                    description='Map from UniProt AA to the new domain positions',
                    epilog='gurdeep')

parser.add_argument('i', help='Input gzipped HMMSearch output file')
parser.add_argument('a', help='Input alignment file')
parser.add_argument('d', help='Name of domain to search for')
parser.add_argument('s', help='Start of domain position')
parser.add_argument('e', help='End of domain position')
args = parser.parse_args()
INPUT_FILE = args.i
ALN_FILE = args.a
DOMAIN = args.d
START_DOMAIN = int(args.s)
END_DOMAIN = int(args.e)

# Map domain position to alignment position
domain2aln = {}
for line in open('../pfam/'+DOMAIN+'.hmm', 'r'):
	if line.split()[-3:] == ['-', '-', '-']:
		domain2aln[int(line.split()[0])] = int(line.split()[-5])

# print (domain2aln)
# sys.exit()

# Map UniProt position to alignment position
acc2aln = {}
for line in open(ALN_FILE, 'r'):
	if line[0] == '>':
		acc = line.split('|')[1]
		if acc not in acc2aln: acc2aln[acc] = {}
		print (line)
		start, end = line.split('|')[3].split('-')
		if start == 'start':
			start = line.split('|')[4].split()[0]
			start = int(start)
		else:
			start = int(start)
		continue
	counter = 0
	for aln_pos in range(0, len(line.rstrip())):
		if line[aln_pos] not in ['-', '.']:
			acc2aln[acc][start+counter] = aln_pos + 1
			counter += 1

# print (acc2aln['Q15418'])
# sys.exit()
flag = 0
dic_kinase = {}
# for line in gzip.open('humanKinasesHitsSplitHmmsearchTrimmed.txt.gz', 'rt'):
for line in gzip.open(INPUT_FILE, 'rt'):
	if len(line.split()) == 0:
		continue
	if 'Domain annotation for each sequence (and alignments):' in line:
		flag = 1
		continue
	if flag == 0:
		continue
	if line[:2] == '>>':
		acc = line.split('>>')[1].lstrip().rstrip().rstrip('\n')
		acc = acc.split()[0]
		if acc not in dic_kinase: dic_kinase[acc] = ''
		# if len(dic_kinase) == 2: break
		# print (acc)
		continue
	if line.split()[0] == DOMAIN:
		start_pfam = int(line.split()[1])
		seq_pfam = line.split()[2]
		end_pfam = int(line.split()[3].rstrip('\n'))
		continue
	if line.split()[0] == acc:
		start_acc = line.split()[1]
		if start_acc == '-':
			continue
		start_acc = int(start_acc)
		seq_acc = line.split()[2]
		end_acc = int(line.split()[3].rstrip('\n'))
		for pfam_aa, acc_aa in zip(seq_pfam, seq_acc):
			if pfam_aa not in ['.', '-'] and acc_aa not in ['.', '-']:
				print (acc, start_acc, start_pfam, pfam_aa, acc_aa)
				if start_pfam >=START_DOMAIN and start_pfam <= END_DOMAIN:
					dic_kinase[acc] += acc + '\t'
					dic_kinase[acc] += str(acc_aa) + '\t' + str(start_acc) + '\t'
					dic_kinase[acc] += str(pfam_aa) + '\t' + str(start_pfam) + '\t'
					# dic_kinase[acc] += str(domain2aln[start_pfam]) + '\n'
					if acc.split('|')[1] in acc2aln:
						if start_acc in acc2aln[acc.split('|')[1]]:
							dic_kinase[acc] += str(acc2aln[acc.split('|')[1]][start_acc]) + '\n'
						else:
							dic_kinase[acc] += '-\n'
					else:
						dic_kinase[acc] += '-\n'
				start_acc += 1
				start_pfam += 1
			elif pfam_aa not in ['.', '-']:
				start_pfam += 1
			elif acc_aa not in ['.', '-']:
				start_acc += 1
		# print (dic_kinase)
		continue
# print (text)
text = '#Acc\tUniProt_AA\tUniProt_Pos\tPfam_AA\tPfam_Pos\n'
for acc in dic_kinase:
	# print (dic_kinase[acc])
	text += dic_kinase[acc]

gzip.open(INPUT_FILE.split('.')[0] + 'Mappings.tsv.gz', 'wt').write(text)