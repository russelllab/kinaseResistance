from Bio import SeqIO

records = SeqIO.parse("muscleAlignment2.fasta", "fasta")
count = SeqIO.write(records, "muscleAlignment2.clustal", "clustal")
print("Converted %i records" % count)
