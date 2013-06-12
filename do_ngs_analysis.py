#!/usr/bin/env python


# prepare raw NGS data for analysis
from pyphylogenomics import NGS;

ionfile = "ionrun.fastq";
#NGS.prepare_data(ionfile);


# find reads matching target genes using BLAST
from pyphylogenomics import BLAST;

query_seqs = "data/modified/wrk_ionfile.fasta";
genome = "target_genes.fasta";
#BLAST.blastn(query_seqs, genome); 

blast_table = "data/modified/wrk_ionfile_blastn_out.csv";
sbj_db      = "data/modified/wrk_ionfile.fastq";
NGS.parse_blast_results(blast_table, sbj_db);

