'''
=====
BLAST
=====

Prepares data and executes BLAST commands to create BLAST database.
Performs blastn of user sequences against genomic sequences.
'''

import re;
import subprocess;
from Bio import SeqIO;
from Bio.SeqRecord import SeqRecord;
from Bio.Blast.Applications import NcbiblastnCommandline;



def get_cds(genes, cds_file):
	'''
	Writes a FASTA file containing CDS sequences: ``pulled_seqs.fa``
	
	``genes`` argument is a list of gene IDs.
	``cds_file`` is the file name of predicted CDS for a species in FASTA format.

	'''

	records = []; # To store the sequences that matched.
	for seq_record in SeqIO.parse(cds_file, "fasta"):
	    this_id = re.sub("-TA$", "", seq_record.id);
	    if this_id in genes:
	        records.append(SeqRecord(seq_record.seq, id=this_id));

	SeqIO.write(records, open("pulled_seqs.fa", "w"), "fasta") 
	print len(records), " sequences were written to file pulled_seqs.fa";



def makeblastdb(genome, mask=False):
	'''
	Creates a BLAST database from a genome in FASTA format and
	optionally eliminates low-complexity regions from the sequences.
	'''

	if mask == True:
	    command = 'dustmasker -in '+ genome + ' -infmt fasta -parse_seqids '
	    command += '-outfmt maskinfo_asn1_bin -out ' + genome + '_dust.asnb'
	    print "masking low_complexity regions..."
	    p = subprocess.check_output(command, shell=True) # identifying low-complexity regions.
	    print p
	    
	    command = 'makeblastdb -in ' + genome + ' -input_type fasta -dbtype nucl '
	    command += '-parse_seqids -mask_data ' + genome + '_dust.asnb '
	    command += '-out ' + genome + ' -title "Whole Genome without low-complexity regions"'
	    print "creating database..."
	    p = subprocess.check_output(command, shell=True) # Overwriting the genome file.
	    print p

	else:
	    command  = 'makeblastdb -in ' + genome + ' -input_type fasta -dbtype nucl '
	    command += '-parse_seqids -out ' + genome + ' -title "Whole Genome unmasked"'
	    print "creating database..."
	    p = subprocess.check_output(command, shell=True) 
	    print p



def blastn(query_seqs, genome):
	'''
	Performs a BLASTn of user's sequences against a genome. It will create a
	BLAST database from the genome file first.

	``seqs`` argument is a FASTA file of users sequences.
	``genome`` argument is a FASTA file of a species genome.
	'''

	# do a BLAST database first
	makeblastdb(genome, mask=False);
	blast_out = query_seqs.split(".")[0]+"_blastn_out.csv" # Name of the output file.

	output  = "blasting: " + query_seqs.split("\\")[-1] + " "
	output += "against db = " + genome.split("\\")[-1]
	print output;

	cline1 = NcbiblastnCommandline(query=query_seqs, db=genome, evalue=0.00001, strand="both",
	                                out=blast_out, num_threads=2, outfmt=10)
									#db_soft_mask=11, out=blast_out, num_threads=2, outfmt=10)
	print cline1
	stdout, stderr = cline1()
	print "\nBLASTn finished!"
	print "The BLAST results were written in to the file ", blast_out

