from Uploader import save_seq;
from DB import create_database;

from Bio import SeqIO;


host = "localhost";
user = "root";
passwd = "mysqlboriska";
db = "pygenomics";

conn = create_database(host, user, passwd, db);

f = open("fasta.list", "r");
lines = f.readlines();

for line in lines:
	line = line.strip();

	for seq_record in SeqIO.parse(line, "fasta"):
		if "NSG" not in str(seq_record.id):
			tmp = str(seq_record.id).split("_");

			geneCode = tmp[0];
			geneName = "NULL";
			code = tmp[1];
			sequence = str(seq_record.seq);

			# save all sequences into goodGenes table
			save_seq(geneCode, geneName, code, sequence, "goodGenes", conn);

