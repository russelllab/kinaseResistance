#!/usr/bin/env python3.10
# coding: utf-8

"""
This script takes alingnment as input
and returns a new alignment without gaps
"""

import argparse
import sys
from Bio import AlignIO
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.AlignIO import ClustalIO


def remove_empty_cols(alignment_file):
    '''
    Removes empty columns from alignment
    '''
    dic_alignment = {}
    for line in open(alignment_file, "r"):
        if line.startswith("CLUSTAL") or line.startswith(" ") \
        or line.startswith("\n") or line.startswith("//"):
            continue
        else:
            # print (line.split())
            name = line.split()[0].strip()
            seq = line.split()[1].strip()
            if name not in dic_alignment:
                dic_alignment[name] = []
            dic_alignment[name] += list(seq)
    
    gap_cols = []
    # print (dic_alignment[name])
    for i in range(len(dic_alignment[name])):
        gap = 0
        for name, seq in dic_alignment.items():
            if seq[i] == "-":
                gap += 1
        if gap == len(dic_alignment):
            gap_cols.append(i)

    print (f"Number of columns to be removed: {len(gap_cols)}")
    fasta_alignment = ''
    for name, seq in dic_alignment.items():
        fasta_alignment += f">{name}\n"
        for i in range(len(seq)):
            if i not in gap_cols:
                fasta_alignment += seq[i]
        fasta_alignment += "\n"
    
    fasta_file = alignment_file.split(".")[0] + ".fasta"
    open(fasta_file, "w").write(fasta_alignment)

    # Read the FASTA alignment
    # fasta_alignment = AlignIO.read(fasta_file, "fasta")
    fasta_alignment = MultipleSeqAlignment(SeqIO.parse(fasta_file, "fasta"))
    records = list(SeqIO.parse(fasta_file, "fasta"))

    # Path to output Clustal file
    clustal_file = alignment_file.split(".")[0] + "_no_gaps.aln"

    # Calculate the maximum sequence name length
    max_name_length = max(len(record.id) for record in records)

    # Write the Clustal alignment
    with open(clustal_file, "w") as output_handle:
        for record in records:
            # Write the sequence name
            output_handle.write("{:<{width}} {}\n".format(record.id, record.seq, width=max_name_length))

    # Write the Clustal alignment
    # with open(clustal_file, "w") as output_handle:
    #     SeqIO.write(fasta_alignment, output_handle, "clustal")

    # Write the Clustal alignment
    # AlignIO.write(fasta_alignment, clustal_file, "clustal")

def checkCLUSTAL(alignment):
    '''
    Checks if alignment file is in CLUSTAL format
    '''
    with open(alignment, "r") as f:
        first_line = f.readline()
        return first_line.startswith("CLUSTAL")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("alignment", help="alignment file in CLUSTAL format")
    args = parser.parse_args()
    alignment = args.alignment

    if checkCLUSTAL(alignment) == False:
        sys.exit("Error: alignment file is not in CLUSTAL format")
    remove_empty_cols(alignment)
    