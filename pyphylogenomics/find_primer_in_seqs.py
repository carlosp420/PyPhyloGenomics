from Bio import SeqIO;
from Bio.Seq import Seq;
from Bio.Alphabet import IUPAC;

import re;
import sys;
import string;


## it assumes that there is a barcode sequence at the start of the read
## if there is no barcode at the start then use 'bp numer of barcodes' = 0
## it will look for primers at the start of the read only

if len(sys.argv) < 4:
	print "Error, need to enter 'primer list file'  'fastq file' and 'barcode bp number' as arguments\n"; 
	print "command.py primerlist.fas fastafiletoprocess.fas 'bp number of barcodes'\n";
	sys.exit();


# get primers 
primer_list = sys.argv[1].strip();

fastaFile = sys.argv[2].strip();
fastaFile_root = re.sub("\.f.+$", "", fastaFile, flags=re.IGNORECASE);

# number of base pairs of the barcode
barcode_bp = sys.argv[3].strip(); 



# -------------------------
# @arguments: degenerated primer;  sequence read
# @brief: iterate through primer's degenerated IUPAC and try to find it in sequence read
# @return TRUE on success
def find_primer_in_seq(primer, seq, barcode_bp):
	found = "false";

	read_seq = str(seq.seq);
	read_seq = read_seq.upper();

	while (found != "true" ):
		primer_seq = str(primer.seq)
		primer_seq = primer_seq.upper();

		# use degenerated code, 
		primer_seq = str(primer_seq).replace('R', '[AG]');
		primer_seq = str(primer_seq).replace('Y', '[CT]');
		primer_seq = str(primer_seq).replace('S', '[GC]');
		primer_seq = str(primer_seq).replace('W', '[AT]');
		primer_seq = str(primer_seq).replace('K', '[GT]');
		primer_seq = str(primer_seq).replace('M', '[AC]');
		primer_seq = str(primer_seq).replace('B', '[CGT]');
		primer_seq = str(primer_seq).replace('D', '[AGT]');
		primer_seq = str(primer_seq).replace('H', '[ACT]');
		primer_seq = str(primer_seq).replace('V', '[ACG]');
		primer_seq = str(primer_seq).replace('N', '[ACGT]');

		my_regex = '^.{' + str(barcode_bp) + '}' + primer_seq;
		m = re.search(my_regex, read_seq);

		if m is not None:
			found = "true";
			return "TRUE";

		# reversecomplement
		primer_seq = str(primer.seq.reverse_complement())
		primer_seq = primer_seq.upper();

		primer_seq = str(primer_seq).replace('R', '[AG]');
		primer_seq = str(primer_seq).replace('Y', '[CT]');
		primer_seq = str(primer_seq).replace('S', '[GC]');
		primer_seq = str(primer_seq).replace('W', '[AT]');
		primer_seq = str(primer_seq).replace('K', '[GT]');
		primer_seq = str(primer_seq).replace('M', '[AC]');
		primer_seq = str(primer_seq).replace('B', '[CGT]');
		primer_seq = str(primer_seq).replace('D', '[AGT]');
		primer_seq = str(primer_seq).replace('H', '[ACT]');
		primer_seq = str(primer_seq).replace('V', '[ACG]');
		primer_seq = str(primer_seq).replace('N', '[ACGT]');

		my_regex = '^.{' + str(barcode_bp) + '}' + primer_seq;
		m = re.search(my_regex, read_seq);
		if m is not None:
			found = "true";
			return "TRUE";


		return "FALSE";

# -------------------------
def main():
	for primer in SeqIO.parse(primer_list, "fasta"):
		# foreach primer try to find them in fasta file
		for fastaFile_seq in SeqIO.parse(fastaFile, "fastq"):
			print "Doing ", fastaFile_seq.id;
			find_degenerated = "";
			find_degenerated = find_primer_in_seq(primer, fastaFile_seq, barcode_bp);
			if find_degenerated == "TRUE":
				print "Match ", primer.id, fastaFile_seq.id, len(str(fastaFile_seq.seq));
				print ">" + fastaFile_seq.id + "_matched_" + primer.id + "\n" + fastaFile_seq.seq;

				# write fastq sequence to file
				out_handle = open(fastaFile_root + "_" + re.sub("_[F|R]", "", primer.id) + ".fastq", "a");
				SeqIO.write(fastaFile_seq, out_handle, "fastq");
				out_handle.close();
		
		

main();
