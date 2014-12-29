#!/usr/bin/env python

'''
=======
OrthoDB
=======

This module uses the table ``OrthoDB6_Arthropoda_tabtext.csv``
downloaded from OrthoDB (the database of orthologous groups)
ftp://cegg.unige.ch/OrthoDB6 to extract certain features and objects
described below.
'''

import argparse


class OrthoDB:
    def __init__(self, in_file, species_name):
        self.in_file = in_file
        self.species_name = species_name
        self.single_copy_genes = self._get_single_copy_genes()
        self.genes = self._get_single_copy_genes()

    def _get_single_copy_genes(self):
        '''
        Returns a list of single-copy genes for a species given by user.

        The species name should be in the same format as stated in the input
        file *in_file*.
        '''
        dictio = self._copies_per_gene()

        ids = dict()
        for key in dictio:
            if self.species_name in key and dictio[key] == 1:
                ortho_id = key[2]
                gene = key[1]

                if ortho_id not in ids:
                    ids[ortho_id] = [gene]
                else:
                    ids[ortho_id].append(gene)
            else:
                pass

        genes = list()
        for k, v in ids.items(): # for k, v in ids.iteritems()
            # we want only single copy genes
            if len(v) < 2:
                genes.append(v[0])
        return(genes)

    def _copies_per_gene(self):
        '''
        *Internal function*

        Creates a dictionary with the number of copies per gene and per species.

        The dictionary has as keys tuples=(species_name, gene_name) and
        values=integer that represent the number of copies of that gene
        in the species.
        '''
        dictio = {}

        handle = open(self.in_file, "rb")
        handle.next()  # skip header
        for line in handle:
            line = line.split('\t')
            specie = line[4]
            gene = line[3]
            o_id = line[1]
            if (specie, gene, o_id) not in dictio:
                dictio[(specie, gene, o_id)] = 1
            else:
                dictio[(specie, gene, o_id)] += 1

        handle.close()
        return(dictio)

    def _single_copy_in_species(self, gene_name):
        '''
        Returns the species where the gene_name given by user is a single-copy gene.

        The gene name should be in the same format as stated in the input file
        *in_file*.
        '''
        print("\nLooking for species with single-copy gene: " + str(gene_name))

        dictio = self._copies_per_gene()

        species = list()
        for key in dictio:
            if gene_name in key and dictio[key] == 1:
                print("Found species: " + str(key[0]))
                species.append(str(key[0]))
            else:
                pass

        print("Found " + str(len(species)) + " species.")
        self.species = species
        return(species)

    def _copies_per_gene_table(self, out_file):
        '''
        Stores the number of copies a gene has per species in a file.

        The output file is a table, where each row is a gene, and each columm
        is a species. The values in the table are the number of gene copies.
        Note: the first row in the otput table is the number of single-copy
        genes in a given species.
        '''
        print("Parsing input file...")

        dictio = self._copies_per_gene()

        genes = list(self.genes)
        genes.sort()
        species = list(self.species)
        species.sort()

        # Writing header of output file.
        out_file_handle = open(out_file, "w")
        out_file_handle.write("Genes\t" + "\t".join(species) + '\n')

        nr = []  # To store number of single-copy genes per specie.
        for sp in species:
            nr.append(0)
            for gn in genes:
                if (sp, gn) in dictio and dictio[(sp, gn)] == 1:
                    nr[-1] += 1  # This is a single-copy gene. Count it!
                else:
                    pass

        nr = [str(i) for i in nr]
        out_file_handle.write("Nr_single_copy_genes\t" + "\t".join(nr) + '\n')

        for gn in genes:
            row = [gn]
            for sp in species:
                if (sp, gn) in dictio:
                    row.append(str(dictio[(sp, gn)]))
                else:
                    row.append(str(0))

            out_file_handle.write("\t".join(row) + '\n')

        out_file_handle.close()

        print("\nOUTPUT FILE WAS GENERATED!")
