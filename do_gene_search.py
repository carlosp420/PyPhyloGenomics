#!/usr/bin/env python

import os;
from pyphylogenomics import OrthoDB
from pyphylogenomics import BLAST


"""
We will find all single-copy genes for the silk moth Bombyx mori using the table
from OrthoDB as input file:
"""
in_file = 'data/OrthoDB6_Arthropoda_tabtext.csv'
genes = OrthoDB.single_copy_genes(in_file, 'Bombyx mori')


"""
Pull all sequences for our gene IDs from the CDS file and write them to a file
pulled_seqs.fa:
"""
cds_file = "data/silkcds.fa"
BLAST.get_cds(genes, cds_file)
print "File moved to data/pulled_seqs.fa"
os.rename("pulled_seqs.fa", "data/pulled_seqs.fa")


"""
Do a BLASTn of the sequences against the Bombyx mori genome. The input arguments
are your file containing the sequences for single-copy genes (pulled_seqs.fa) 
and your file with the genome of Bombyx mori which is in FASTA format (silkgenome.fa).
"""
BLAST.blastn('data/pulled_seqs.fa', 'data/silkgenome.fa')
print "File moved to data/pulled_seqs_blastn_out.csv"
os.rename("pulled_seqs_blastn_out.csv", "data/pulled_seqs_blastn_out.csv")
