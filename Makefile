# This is the procedure to reproduce the analyses in the PyPhyloGenomics manuscript
#

## Get data files
# Orthodb file
data/OrthoDB6_Arthropoda_tabtext.csv: 
	wget ftp://cegg.unige.ch/OrthoDB6/OrthoDB6_Arthropoda_tabtext.gz
	mkdir -p data
	mv OrthoDB6_Arthropoda_tabtext.gz data/.
	gunzip data/OrthoDB6_Arthropoda_tabtext.gz
	mv data/OrthoDB6_Arthropoda_tabtext data/OrthoDB6_Arthropoda_tabtext.csv

# Bombyx mori CDS
data/silkcds.fa:
	wget ftp://silkdb.org/pub/current/Gene/Glean_genes/silkworm_glean_cds.fa.tar.gz
	mv silkworm_glean_cds.fa.tar.gz data/.
	tar -zxvf data/silkworm_glean_cds.fa.tar.gz -C data/


# Bombyx mori genome 
data/silkworm_genome_v2.0.fa.tar.gz: 
	curl -0 ftp://silkdb.org/pub/current/Genome/silkworm_genome_v2.0.fa.tar.gz > data/silkworm_genome_v2.0.fa.tar.gz

data/silkgenome.fa: data/silkworm_genome_v2.0.fa.tar.gz
	tar -zxvf data/silkworm_genome_v2.0.fa.tar.gz -C data/


# Danaus genome
data/Dp_genome_v2.fasta:
	# The original link is not working!
	# curl -0 http://monarchbase.umassmed.edu/download/Dp_genome_v2.fasta.gz > data/Dp_genome_v2.fasta.gz
	# Get a copy from own server
	curl -0 http://nymphalidae.utu.fi/cpena/etc/Dp_genome_v2.fasta.gz > data/Dp_genome_v2.fasta.gz
	gunzip data/Dp_genome_v2.fasta.gz

# Heliconius genome
data/Heliconius_genome.fa:
	curl -0 ftp://ftp.ensemblgenomes.org/pub/metazoa/release-18/fasta/heliconius_melpomene/dna/Heliconius_melpomene.Hmel1.18.dna.toplevel.fa.gz > data/Heliconius_genome.fa

gene_search: data/OrthoDB6_Arthropoda_tabtext.csv data/silkcds.fa data/silkgenome.fa data/Dp_genome_v2.fasta data/Heliconius_genome.fa
	echo "Done"

