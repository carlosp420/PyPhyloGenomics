import sys;
from pyphylogenomics import NGS;

fastq_file   = "output/index_IonADA_656_gene_cad_1.fastq";
index_length = 8;
min_quality  = 20; # optional
percentage   = 70; # optional
min_length   = 50; # optional
NGS.assembly(fastq_file, index_length, min_quality, percentage, min_length);

