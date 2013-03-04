from Bio import SeqIO;
from Bio.Seq import Seq;
from Bio.Alphabet import IUPAC;

import re;
import sys;
import string;
import os.path


if len(sys.argv) < 3:
	print "\nError, need to enter barcodes list file and fasta file as arguments"; 
	print "Allows 1 mismatch. Our barcodes have a difference of 2 base pairs.\n";
	print "command.py barcodes_.fas fastqFiletoProcess.fastq\n";
	sys.exit();

# get primers 
barcodes_list = sys.argv[1].strip();

fastqFile = sys.argv[2].strip();
fastqFile_root = re.sub(".fastq", "", fastqFile);

fastq_id_list = [];


# -------------------------
def levenshtein(a,b):
    "Calculates the Levenshtein distance between a and b."
    n, m = len(a), len(b)
    if n > m:
        # Make sure n <= m, to use O(min(n,m)) space
        a,b = b,a
        n,m = m,n
        
    current = range(n+1)
    for i in range(1,m+1):
        previous, current = current, [i]+[0]*n
        for j in range(1,n+1):
            add, delete = previous[j]+1, current[j-1]+1
            change = previous[j-1]
            if a[j-1] != b[i-1]:
                change = change + 1
            current[j] = min(add, delete, change)
            
    return current[n]

# -------------------------
# @arguments:  barcodes;  sequence read
# @brief: iterate through primer's degenerated IUPAC and try to find it in sequence read
# @return TRUE on success
def find_barcode_in_seq(barcode, seq):
	found = "false";

	while (found != "true" ):
		barcode_seq = str(barcode.seq) 
		barcode_read = str(seq.seq)[0:8];

		#compare levenshtein
		distance = levenshtein(barcode_seq.upper(), barcode_read.upper());

		# if distance < 2 :
		if distance < 2:
			# accept
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
	filename = "reads_with_no_barcodes.fastq";

	# avoid overwritting output file
	while os.path.exists(filename) == True:
		count = re.sub("reads_with_no_barcodes","", filename)
		count = re.sub(".fastq","", count)
		if count == "":
			count = 1;
		else:
			count = int(count) + 1;
		tmp = str(count) + ".fastq";
		filename = re.sub(".fastq", tmp, filename);

	out_handle = open(filename, "w");
	for fastqFile_seq in SeqIO.parse(fastqFile, "fastq"):
		if fastqFile_seq.id not in fastq_id_list:
			SeqIO.write(fastqFile_seq, out_handle, "fastq");
	out_handle.close();


# -------------------------
def main():

	for seq_record in SeqIO.parse(barcodes_list, "fasta"):
		# foreach barcode try to find them in fastq file
		for fastqFile_seq in SeqIO.parse(fastqFile, "fastq"):
#print "Doing ", fastqFile_seq.id;
			find_barcode = "";
			find_barcode = find_barcode_in_seq(seq_record, fastqFile_seq);
			if find_barcode == "TRUE":
				print "##########Match ", seq_record.id, fastqFile_seq.id, len(str(fastqFile_seq.seq));

				output_handle = open(str(seq_record.id) + fastqFile_root + ".fastq", "a");
				SeqIO.write(fastqFile_seq, output_handle, "fastq");
				output_handle.close();
	
		print fastq_id_list;

	pull_out_reads(fastq_id_list);
	"""
	pull out all reads that are not in this list
	"""
		
		

main();
