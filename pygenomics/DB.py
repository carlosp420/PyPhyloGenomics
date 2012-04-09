# This code is part of the PyGenomics distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""

Collection of modules for dealing with genomic data in Python.

Interact with MySQL database
"""

import MySQLdb;


"""
creates database "pygenomics" if does not exist
"""
def create_database(host, user, passwd, db):
	try:
		mysql = MySQLdb.connect(unix_socket="/tmp/mysql.sock", host=host, user=user, passwd=passwd, db=db);
	except MySQLdb.Error:
		mysql = MySQLdb.connect(unix_socket="/tmp/mysql.sock", host=host, user=user, passwd=passwd);
		mysql.query("create database pygenomics");
		mysql = MySQLdb.connect(unix_socket="/tmp/mysql.sock", host=host, user=user, passwd=passwd, db=db);

	mysql.query("set names utf8");


	"""create table "good_genes" to host all sequences for genes that have to
		be sequenced because are < 300bp"""

	table = "CREATE TABLE IF NOT EXISTS `goodGenes` ( \
				`id` smallint(5) unsigned NOT NULL AUTO_INCREMENT, \
				`geneCode` varchar(255) DEFAULT NULL, \
				`geneName` varchar(255) DEFAULT NULL, \
				`code` varchar(255) DEFAULT NULL, \
				`sequence` text DEFAULT NULL, \
				`timestamp` datetime NOT NULL DEFAULT '0000-00-00 00:00:00', \
				PRIMARY KEY (`id`), \
				UNIQUE KEY `id` (`id`) \
				) ENGINE=MyISAM AUTO_INCREMENT=1 DEFAULT CHARSET=utf8"; 
	mysql.query(table);

	return mysql;
