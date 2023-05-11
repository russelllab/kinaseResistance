#!/usr/bin/env python3

'''
Convert FASTA to CLUSTAL
'''

from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(
                    prog='fasta2clustal',
                    description='Converts aln in FASTA format to CLUSTAL format',
                    epilog='gurdeep')

parser.add_argument('i', help='Input alignment in FASTA format')
parser.add_argument('o', help='Output alignment in CLUSTAL format')
args = parser.parse_args()
input_file = args.i
output_file = args.o

records = SeqIO.parse(input_file, "fasta")
count = SeqIO.write(records, output_file, "clustal")
print("Converted %i records" % count)
