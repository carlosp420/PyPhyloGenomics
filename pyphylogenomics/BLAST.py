'''
=====
BLAST
=====

Prepares data and executes BLAST commands to create BLAST database

@input:  FASTA file containing genomic sequences in working directory
@param:  FASTA filename 
@output: BLAST database in working directory
'''

import re;
import subprocess;
from Bio import SeqIO;



def get_cds(genes, cds_file):
	'''
	Writes a FASTA file containing CDS sequences: ``pulled_seqs.fa``
	
	``genes`` argument is a list of gene IDs.
	``cds_file`` is the file name of predicted CDS for a species in FASTA format.

	'''
	outfile = open("pulled_seqs.fa", "w");

	i = 0;
	for seq_record in SeqIO.parse(cds_file, "fasta"):
		this_id = re.sub("-TA$", "", seq_record.id);
		if this_id in genes:
			outfile.write(">" + this_id + "\n");
			outfile.write(str(seq_record.seq) + "\n");
			i = i + 1;
	outfile.close();

	print i, " sequences were written to file pulled_seqs.fa";


def makeblastdb(FASTA_file):
	'''
	Creates a BLAST database from a genome in FASTA format.
	'''

	command  = 'makeblastdb -in ' + FASTA_file;
	command += ' -dbtype nucl -parse_seqids -input_type fasta'; 
	p = subprocess.check_output(command, shell=True) 
	print p
	
	command  = 'blastdb_aliastool -dblist "' + FASTA_file + '" ';
	command += '-dbtype nucl -out ' +  FASTA_file + ' -title "' + FASTA_file + '"';
	p = subprocess.check_output(command, shell=True)
	print p

	command  = 'makembindex -input ' + FASTA_file;
	command += ' -iformat fasta -output ' + FASTA_file;
	p = subprocess.check_output(command, shell=True)
	print p


def blastn(seqs, genome):
	'''
	Performs a BLASTn of user's sequences against a genome. It will create a
	BLAST database from the genome file first.

	``seqs`` argument is a FASTA file of users sequences.
	``genome`` argument is a FASTA file of a species genome.
	'''

	# do a BLAST database first
	makeblastdb(genome);

	# ok, now do the BLAST
	# @ input is the BLAST database created by the command ``makeblastdb(genome)``
	# @ output is the output file as a table in a csv file.
	# TODO:
		# This was written by Chris, but it should be corrected/changed.
		# Make sure that the BLAST works and that the table is in the right format
		# Make sure that it works ok, and returns useful messages to the user (prints to screen)

    print "blasting: ", blast_query
    print "against db = ",blast_db
    cline1 = NcbiblastxCommandline(cmd=blast_exe, query=blast_query, db=blast_db, evalue=0.00001, outfmt=6, out=blast_out, num_threads=2)
    return_code1 = subprocess.call(str(cline1), shell=(sys.platform!="linux")) # this is calling the blast
    print return_code1
    print "blast finished"


