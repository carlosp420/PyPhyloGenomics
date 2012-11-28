"""
MODULE
------------------------------------------------------------------------------
execute BLAST commands to create BLAST database
------------------------------------------------------------------------------
@input:  FASTA file containing genomic sequences in working directory
@param:  FASTA filename 
@output: BLAST database in working directory
"""

# Python recommends to use now subprocess.Popen() instead of os.popen()
import subprocess


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
