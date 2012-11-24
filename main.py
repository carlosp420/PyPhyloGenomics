from Uploader import save_seqs;
from DB import create_database;
from DB import make_blastdb;

from Bio import SeqIO;


host = "localhost";
user = "root";
passwd = "mysqlboriska";
db = "pygenomics";

conn = create_database(host, user, passwd, db);

def save_my_seqs():
	f = open("fasta.list", "r");
	lines = f.readlines();

	for line in lines:
		line = line.strip();

		for seq_record in SeqIO.parse(line, "fasta"):
			if "NSG" not in str(seq_record.id):
				tmp = str(seq_record.id).split("_");

				geneName = tmp[0];
				geneCode = tmp[1];
				code = tmp[2];
				sequence = str(seq_record.seq);
				sequence = sequence.strip("?");

				# save all sequences into goodGenes table
				save_seqs(geneCode, geneName, code, sequence, "goodGenes", conn);

#save_my_seqs();
make_blastdb(conn);
