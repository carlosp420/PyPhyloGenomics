import sys;
from Bio import SeqIO;

input_file = sys.argv[1];

out = sys.argv[2];
out_handle = open(out, "w");

for seq_in in SeqIO.parse(input_file, "fastq"):
	print seq_in.id;
	SeqIO.write(seq_in, out_handle, "fasta");

out_handle.close()
