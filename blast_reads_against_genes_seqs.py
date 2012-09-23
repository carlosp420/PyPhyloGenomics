#!/usr/bin/python

from Bio import SeqIO;
import sys;
import os;
import time;
import string;
import subprocess, shlex;
import re;
from Bio.Blast.Applications import NcbiblastnCommandline

#############
if len(sys.argv) < 2:
	print "You need to enter as arguments: (i) your 'reads_file.fas' and (ii) your database name 'all_genes.fas'";
	print "command.py reads_file.fas all_genes.fas";
	# Input file is a fasta file including several DNA sequences
	sys.exit();

exon_length = int("40");
#exon_length = "300";
filename = sys.argv[1];
filename_root = re.sub("\.f.{1,4}$", "", filename);

db_name = sys.argv[2];
db_name_root = db_name.replace(".fas", "");
fileoutput = open("reads_against_seqs_blast_out_" + filename_root + ".csv", "w");


def extract_table_hit(result):
	result = result.split("\n");
	for line in result:
		line = line + "\n";
		return line;
		break;

def pull_sequence(result):
	result = result.split("\t");
	for seq_record in SeqIO.parse(db_name, "fasta"):
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
	#print result[3];
	if int(result[3]) > exon_length:
		return result[3];

def check_if_need_reverse_complement(result):
	result = result.split("\t");
	if int(result[8]) > int(result[9]):
		# subject start > subject end
		# need to be reverse complemented
		return "reverse";

fileoutput.write("Query\tSubject\t%id\tAlignment length\tmismatches\tgaps\tq.start\tq.end\ts.start\ts.end\te.value\tbitscore\n");


for seq_record in SeqIO.parse(filename, "fastq"):
	print "doing... ", seq_record.id;
	if os.path.exists("query.fas"):
		os.unlink("query.fas");

	out = open("query.fas", "w");
	string = ">" + str(seq_record.id) + "\n" + str(seq_record.seq);
	#print string;
	out.write(string);
	out.close();


#command = 'blastall -p blastn -i query.fas -d "Dp_genome_v2.fasta" -e 0.01 -m 8';
	command = 'blastn -query query.fas -db ' + db_name + ' -task blastn -dust no -outfmt 6 -num_alignments 4 -num_descriptions 4 -index_name ' + db_name_root + '.00.idx';
	p = os.popen(command);

	result = p.read();
	result = extract_table_hit(result);
	if len(result) > 1:
		# check for alignments > 300bp
		length_exon = int("1");
		length_exon = check_length_exon(result);

		if length_exon != None:
			if int(length_exon) > exon_length:
				print "----------";
				print length_exon;
				print "======";
				print result;
				fileoutput.write(result);


				# create file for exon based on primer and barcode
				file_name = filename_root + "_" + re.sub("_B.+$", "", result.split("\t")[1]);
				exon_out = open(file_name + ".fastq", "a");

				"""
				output = ">" + str(seq_record.id) + "\n";
						
				need_to_reverse = "";
				need_to_reverse = check_if_need_reverse_complement(result);

				if need_to_reverse == "reverse":
					output += str(seq_record.seq.reverse_complement()) + "\n";
				else:
					output += str(seq_record.seq) + "\n";

				exon_out.write(output);
				"""

				SeqIO.write(seq_record, exon_out, "fastq");
				exon_out.close();

				"""
				pulled_seq = pull_sequence(result);
				exon_out.write(pulled_seq);
				"""

fileoutput.close()
