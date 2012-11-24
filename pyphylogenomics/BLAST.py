"""
------------------------------------------------------------------------------
execute BLAST commands to create BLAST database
------------------------------------------------------------------------------
"""

import subprocess;
import os;

def makeblastdb():
	command = "makeblastdb -in db.fas -dbtype nucl -parse_seqids -input_type fasta";
	p = os.popen(command);

	command = 'blastdb_aliastool -dblist "db.fas" -dbtype nucl -out db.fas -title "db" ';
	p = os.popen(command);

	command = 'makembindex -input db.fas -iformat fasta -output db';
	p = os.popen(command);
