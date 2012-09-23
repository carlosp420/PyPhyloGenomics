from Bio import SeqIO;
from Bio.Seq import Seq;
from Bio.Alphabet import IUPAC;
import sys;


input_file = sys.argv[1].strip();
output_file = open("output_file", "w");

for seq_record in SeqIO.parse(input_file, "fastq"):
	a = seq_record.reverse_complement();
	a.id = seq_record.id;
	SeqIO.write(a, output_file, "fastq");
	
output_file.close();
