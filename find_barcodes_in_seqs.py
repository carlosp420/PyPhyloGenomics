from Bio import SeqIO;
from Bio.Seq import Seq;
from Bio.Alphabet import IUPAC;

import re;
import sys;
import string;


if len(sys.argv) < 3:
	print "Error, need to enter barcodes list file and fasta file as arguments\n"; 
	print "command.py barcodes_.fas fastafiletoprocess.fas\n";
	sys.exit();

# get primers 
barcodes_list = sys.argv[1].strip();

fastqFile = sys.argv[2].strip();

fastq_id_list = [];

# -------------------------
# @arguments:  barcodes;  sequence read
# @brief: iterate through primer's degenerated IUPAC and try to find it in sequence read
# @return TRUE on success
def find_barcode_in_seq(barcode, seq):
	found = "false";

	while (found != "true" ):
		barcode_seq = str(barcode.seq) 

		m = re.search('^%s' % barcode_seq, str(seq.seq), re.IGNORECASE);
		if m is not None:
			found = "true";

			if str(seq.id) not in fastq_id_list:
				fastq_id_list.append(str(seq.id));
			return "TRUE";

		# reversecomplement
			"""
		barcode_seq = str(barcode.seq.reverse_complement())

		m = re.search(barcode_seq, str(seq));
		if m is not None:
			found = "true";
			return "TRUE";
			"""

		return "FALSE";

# -------------------------
"""
pull out all reads that are not in this list
"""
def pull_out_reads(fastq_id_list):
		for fastqFile_seq in SeqIO.parse(fastqFile, "fastq"):
			if fastqFile_seq.id not in fastq_id_list:
				out_handle = open("reads_with_no_barcodes.fastq", "a");
				SeqIO.write(fastqFile_seq, out_handle, "fastq");
				out_handle.close();


# -------------------------
def main():

	for seq_record in SeqIO.parse(barcodes_list, "fasta"):
		# foreach barcode try to find them in fastq file
		for fastqFile_seq in SeqIO.parse(fastqFile, "fastq"):
#print "Doing ", fastqFile_seq.id;
			find_degenerated = "";
			find_barcode = find_barcode_in_seq(seq_record, fastqFile_seq);
			if find_barcode == "TRUE":
				print "##########Match ", seq_record.id, fastqFile_seq.id, len(str(fastqFile_seq.seq));

				output_handle = open(str(seq_record.id), "a");
				SeqIO.write(fastqFile_seq, output_handle, "fastq");
				output_handle.close();
	
		print fastq_id_list;

	pull_out_reads(fastq_id_list);
	"""
	pull out all reads that are not in this list
	"""
		
		

main();
