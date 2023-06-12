'''
Script to map UniProt seq positions
to Pkinase domain positions
'''
import os, gzip, sys
import argparse, threading
import os, sys, gzip

parser = argparse.ArgumentParser(
                    prog='MapUniProt2Aln',
                    description='Map from UniProt AA to the ALN positions',
                    epilog='gurdeep')

parser.add_argument('a', help='Input alignment file')
parser.add_argument('d', help='Name of domain to search for')
args = parser.parse_args()

ALN_FILE = args.a
DOMAIN = args.d

class KinaseDomain:
	def __init__(self, name, start, end):
		self.name = name
		self.start = start
		self.end = end
		self.aln = ''

	def __str__(self):
		return self.acc + '\t' + str(self.pos) + '\t' + str(self.aln_pos) + '\n'
	
	def name2aln(self, aln2domain):
		counter = 0
		text = ''
		for aln_pos in range(0, len(self.aln.rstrip())):
			if self.aln[aln_pos] not in ['-', '.']:
				uniprotAA = self.aln[aln_pos].upper()
				uniprotPos = self.start+counter
				text += self.name + '\t' + uniprotAA + '\t' + str(uniprotPos) + '\t'
				text += aln2domain[aln_pos+1]['AA'] + '\t' + str(aln2domain[aln_pos+1]['Pos']) + '\t'
				text += str(aln_pos + 1) + '\n'
				counter += 1
		return text


# Map domain position to alignment position
domain2aln = {}
aln2domain = {}
for line in open('../pfam/'+DOMAIN+'.hmm', 'r'):
	if line.split()[-3:] == ['-', '-', '-']:
		domain2aln[int(line.split()[0])] = int(line.split()[-5])
		# aln2domain[int(line.split()[-5])] = int(line.split()[0])
		aln2domain[int(line.split()[-5])] = {'AA': line.split()[-4], 'Pos': int(line.split()[0])}

# print (domain2aln)
# sys.exit()

# Map UniProt position to alignment position
# acc2aln = {}
dic_kinases = {}
for line in open(ALN_FILE, 'r'):
	if line[0] == '>':
		name = line.split('>')[1].split()[0]
		start, end = line.split('|')[3].split('-')
		counter = line.split('|')[4].split()[0]
		if start == 'start':
			start = int(counter)
		else:
			start = int(start) + int(counter) - 1
		if name not in dic_kinases: dic_kinases[name] = KinaseDomain(name, start, end)
		continue
	dic_kinases[name].aln += line.rstrip()

text = '#Name\tUniProt_AA\tUniProt_Pos\tPfam_AA\tPfam_Pos\tAln_Pos\n'
for name in dic_kinases:
		text += dic_kinases[name].name2aln(aln2domain)

gzip.open(DOMAIN + 'Mappings.tsv.gz', 'wt').write(text)