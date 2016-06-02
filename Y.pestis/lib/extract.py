from Bio import SeqIO
import sys

file_in = sys.argv[1]
name = sys.argv[2]
file_out = open(name + "_chr", "w")
for seq in SeqIO.parse(sys.argv[1], "fasta"):
    if name in seq.description:
       	print seq.id
       	SeqIO.write(seq, file_out, "fasta")

