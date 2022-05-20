from Bio import SeqIO

records = SeqIO.parse("muscleAlignment.aln", "fasta")
count = SeqIO.write(records, "muscleAlignment.clustal", "clustal")
print("Converted %i records" % count)
