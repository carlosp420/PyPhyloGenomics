# This code is part of the PyPhyloGenomics distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""

Collection of modules for dealing with genomic data in Python.

"""




__docformat__ ="epytext en" #Don't just use plain text in epydoc API pages!


def save_seqs(geneCode, geneName, code, sequence, table, conn):
	# table needs to be specified by user
	"""
		input: geneCode,
				geneName,
				code,
				sequence,
				timestamp
		actions:
			save to MySQL
	"""
	query = "INSERT INTO " + table + " (geneCode, geneName, code, sequence, timestamp) values ";
	query += "('" + geneCode + "', '" + geneName + "', '" + code + "', '" + sequence + "', ";
	query += "now())";

	print query;
	conn.query(query);



