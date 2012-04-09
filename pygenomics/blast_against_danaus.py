from Bio import SeqIO;
import sys;
import os;
import time;
import subprocess, shlex;
from Bio.Blast.Applications import NcbiblastnCommandline

database_name = "Dp_genome_v2.fasta";

def extract_table_hit(result):
	result = result.split("\n");
	for line in result:
		line = line + "\n";
		return line;
		break;

def pull_sequence(result):
	result = result.split("\t");
	for seq_record in SeqIO.parse(database_name, "fasta"):
		if seq_record.id == result[1]:
			pulled_seq = ">Danaus_plexippus_" + str(seq_record.id) + "\n";
			hit_start = int(result[8]) - 1;
			hit_end   = int(result[9]) + 1;
			print hit_start,  ", " , hit_end 
			if hit_start < hit_end:
				print "forward\n";
				pulled_seq = pulled_seq + str(seq_record.seq[hit_start:hit_end]) + "\n";
			else:
				print "reverse\n";
				seq = seq_record.seq;
				seq = seq[hit_end:hit_start]
				pulled_seq = pulled_seq + str(seq.reverse_complement()) + "\n";

			return pulled_seq;

# function to check for alignments > 300bp
def check_length_exon(result):
	result = result.split("\t");
	if int(result[3]) > 300:
		return result[3];


#############
if len(sys.argv) < 2:
	print "You need to enter input file as argument: \n command.py filename\n";
	# Input file is a fasta file including several DNA sequences
	sys.exit();

filename = sys.argv[1];


fileoutput = open("Bmori_vs_Dplexippus.csv", "w");
fileoutput.write("Query\tSubject\t%id\tAlignment length\tmismatches\tgaps\tq.start\tq.end\ts.start\ts.end\te.value\tbitscore\n");


for seq_record in SeqIO.parse(filename, "fasta"):
	print "doing... ", seq_record.id;
	if os.path.exists("query.fas"):
		os.unlink("query.fas");

	out = open("query.fas", "w");
	out.write(">" + str(seq_record.id) + "\n" + str(seq_record.seq));
	out.close();


#command = 'blastall -p blastn -i query.fas -d "Dp_genome_v2.fasta" -e 0.01 -m 8';
	command = './blastn -query query.fas -db Dp_genome_v2.fasta -task blastn -dust no -outfmt 6 -num_alignments 4 -num_descriptions 4 -index_name Dp_genome_v2.00.idx';
	p = os.popen(command);

	result = p.read();
	result = extract_table_hit(result);
	if len(result) > 1:
		# check for alignments > 300bp
		length_exon = check_length_exon(result);
		if length_exon > 300:
			print result;
			fileoutput.write(result);


			# create file for exon
			exon_out = open(str(seq_record.id) + ".fas", "w");
			exon_out.write(">" + str(seq_record.id) + "\n" + str(seq_record.seq) + "\n");
			pulled_seq = pull_sequence(result);
			exon_out.write(pulled_seq);

fileoutput.close()
