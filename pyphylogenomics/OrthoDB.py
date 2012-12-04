#!/usr/bin/env python

import argparse;
import sys;





def copies_per_gene(in_file):
    '''
	Creates a dictionary with the number of copies per gene and per species.
    
    The dictionary has as keys tuples=(species_name, gene_name) and
    values=integer that represent the number of copies of that gene
    in the species.
    '''


    dictio = {}

    global genes, species
    genes=set()
    species=set()

    in_file = open(in_file, "rb");
    in_file.next() # skip header
    for line in in_file:
        specie = line.strip().split('\t')[4]
        gene = line.strip().split('\t')[3]
        species.add(specie)
        genes.add(gene)
        if (specie,gene) not in dictio:
            dictio[(specie,gene)] = 1
        else:
            dictio[(specie,gene)] += 1

    in_file.close()

    return dictio



def single_copy_genes(in_file, species_name):
    ''' For a species_name given by the user, returns its single-copy genes.

    The species name should be in the same format as stated in the input file
    of the function 'copies_per_gene'.
    '''

    print "\n";
    print "Looking for single-copy genes for species "+ str(species_name) +":";
    
    dictio = copies_per_gene(in_file)
    for key in dictio:
        if species_name in key and dictio[key] == 1:
            print "Found gene: " + str(key[1]);
        else:
            pass

	#in_file.close()
   

def single_copy_in_species(in_file, gene_name):
    ''' For a gene_name given by the user, returns the species where it is
    in a single copy.

    The species name should be in the same format as stated in the input file
    of the function 'copies_per_gene'.
    '''

    print "running..."
    
    dictio = copies_per_gene(in_file)
    for key in dictio:
        if gene_name in key and dictio[key] == 1:
            print key[0]
        else:
            pass

    in_file.close()


def copies_per_gene_table(in_file, out_file):
    ''' 
	This script stores the number of copies a gene has per species in a file.

    The output file is a table, where each row is a gene, and each columm
    is a species. The values in the table are the number of gene copies.
    Note: the first row in the otput table is the number of single-copy
    genes in a given species.
    '''
    
    print "Parsing input file..."

    dictio = copies_per_gene(in_file)
    
    global genes, species
    genes = list(genes)
    genes.sort()
    species = list(species)
    species.sort()

    # Writing header of output file.
    out_file.write("Genes\t"+ "\t".join(species)+'\n')

    nr = [] # To store number of single-copy genes per specie. 
    for sp in species:
        nr.append(0)
        for gn in genes:
            if (sp,gn) in dictio and dictio[(sp,gn)] == 1:
                nr[-1] += 1 # This is a single-copy gene. Count it!
            else:
                pass
    
    nr = [str(i) for i in nr]
    out_file.write("Nr_single_copy_genes\t" + "\t".join(nr) + '\n')

    for gn in genes: 
        row = [gn]
        for sp in species:
            if (sp,gn) in dictio:
                row.append(str(dictio[(sp,gn)]))
            else:
                row.append(str(0))

        out_file.write("\t".join(row) + '\n')


	#in_file.close()
    out_file.close()

    print "\nOUTPUT FILE WAS GENERATED!"        


# -----------------------------------------------------------------------------
def main():
	description = '''This module uses the table OrthoDB6_Arthropoda_tabtext.csv
		downloaded from OrthoDB (the database of orthologous groups) to
		extract certain features and objects described below.  '''

	parser = argparse.ArgumentParser(description=description);
	parser.add_argument('-i', '--input', action='store', 
			metavar='OrthoDB_xyz.csv',
			required=True, dest='in_file',
			help='OrthoDb table as input file');
	parser.add_argument('-o', '--output', action='store', 
			metavar='output_table.csv',
			required=True, dest='out_file',
			help='parsed data as csv file');
	parser.add_argument('-s', '--species', action='store', 
			metavar='Bombyx_mori',
			dest='species_name',
			help='return single-copy genes for given species name (use \
				underscore to join genus and species name)');
	
	args = parser.parse_args();
	
	
	in_file = args.in_file.strip();
	out_file = open(args.out_file, 'w');
	species_name = args.species_name;

	## do default processing
	print "Finding number of copies per gene..."
	copies_per_gene_table(in_file, out_file);

	## if --species argument given
	if species_name != None:
		species_name = species_name.strip();

		# species will be "Bombyx mori"
		species_name = species_name.replace("_", " ");
		single_copy_genes(in_file, species_name);


if __name__ == "__main__":
	main();
