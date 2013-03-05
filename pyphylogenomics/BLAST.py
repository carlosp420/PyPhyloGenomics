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
	Writes a FASTA file containing CDS sequences: ``pulled_seqs.fas``
	
	``genes`` argument is a list of gene IDs.
	``cds_file`` is the file name of predicted CDS for a species in FASTA format.

	'''
	outfile = open("pulled_seqs.fas", "w");

	i = 0;
	for seq_record in SeqIO.parse(cds_file, "fasta"):
		this_id = re.sub("-TA$", "", seq_record.id);
		if this_id in genes:
			outfile.write(">" + this_id + "\n");
			outfile.write(str(seq_record.seq) + "\n");
			i = i + 1;
	outfile.close();

	print i, " sequences were written to file pulled_seqs.fas";


def makeblastdb(FASTA_file):
	command  = 'makeblastdb -in ' + FASTA_file;
	command += ' -dbtype nucl -parse_seqids -input_type fasta'; 
	p = subprocess.Popen(command, shell=True) 
	
	command  = 'blastdb_aliastool -dblist "' + FASTA_file + '" ';
	command += '-dbtype nucl -out ' +  FASTA_file + ' -title "db"';
	p = subprocess.Popen(command, shell=True)

	command  = 'makembindex -input ' + FASTA_file;
	command += ' -iformat fasta -output db'
	p = subprocess.Popen(command, shell=True)
