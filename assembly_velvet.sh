for k in 31 29 27 25 23 21 19 17 15 13 11 9 7 5 3 1
do
	velveth_de test $k -fastq -shortPaired $1
	velvetg_de test  -min_contig_lgth 100
done
