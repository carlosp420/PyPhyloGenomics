from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from operator import itemgetter


def getLargestExon(blast_table, E_value = 0.001, ident = 98, exon_len = 300):
    """
    Returns the largest alignment for each query-sbj match from a blast table. 
    """
    table = open(blast_table, "r") 
    exons = {} # A dictionary, where the key is the (query,sbj) tuple
                # and the value is the largest alignment of all the possible
                # alignments between query and sbj. 

    print "Parsing BLAST table ..."  
    for alignment in table:
        if 'Query' in alignment: continue # in case the input table has header.     

        align_vars = alignment.strip().split(",") # "blast_table" must be a csv file.
        align_vars[3:10] = [int(i) for i in align_vars[3:10]]
        align_vars[3] -= align_vars[5] # replaces the alignment length for the exon length

        if float(align_vars[10]) <= E_value and float(align_vars[2]) >= ident:
            if (align_vars[0],align_vars[1]) in exons:
                exons[(align_vars[0],align_vars[1])][-1] += 1 # Counts the number of alignments with that sbj.
                if align_vars[3] > exons[(align_vars[0],align_vars[1])][3]:
                    exons[(align_vars[0],align_vars[1])][:-1] = align_vars[:]  
            else:
                exons[(align_vars[0],align_vars[1])] = align_vars[:] + [1] # The 1st alignment with that sbj.

    bad = []
    print "Deleting exons below %d nucleotides ..." % exon_len
    for key in exons: # For deleting (query, sbj) alignments whose largest exon is below exon_len threshold. 
        if exons[key][3] < exon_len:
            bad.append(key)
    for key in bad:
        del(exons[key])
         
    print "There are %s exons" % len(exons);
    return exons



def eraseFalsePosi(exons_dict):
    """
    Keeps the query-sbj match with the greatest number of alignments.
    
    From a dictionary generated with the getLargestExon function, where 2 or more
    query-sbj matches shared the same query, eraseFalsePosi keeps the query-sbj
    match with the greatest number of alignments.
    """

    print "Erasing False Positives ..."
    new_exons_dict = {}
    keys = sorted([j + (exons_dict[j][-1],) for j in exons_dict.keys()],
                  key=itemgetter(0,2), reverse=True)
    keys = [j[:2] for j in keys]
    unique = list(set([i[0] for i in keys]))
    for i in unique:
        for key in keys:
            if key[0] == i:
                new_exons_dict.update({key:exons_dict[key]})
                break
        continue

    print "There are %s exons" % len(new_exons_dict);
    return new_exons_dict



def filterByMinDist(genes_loci, MinDist):
    """
    Returns those genes that are separated from its precedent neighbour
    by < MinDist threshold.

    The parameter genes_loci is a list of tuples, where each tuple
    consists of 3 items: a gene_id, start_position and end_position.
    """
    
    genes_loci.sort(key=itemgetter(1,2))
    genes = []
    last_end = genes_loci[0][2]
    for loci in genes_loci[1:]:
        if loci[1] - last_end < MinDist:
            genes.append(loci[0])
            last_end = loci[2]

    return genes
                                              


def wellSeparatedExons(exons_dict, MinDist=810000):
    """
    Keeps the exons whose distance to the following exon is > MinDist.
    
    From a dictionary generated with the getLargestExon function, where 2 or more
    query-sbj matches shared the same sbj, wellSeparatedExons keeps the query-sbj
    match whose distance to the following query-sbj match is greater than MinDist.
    """
    print "Identifying exons separated by %d bases ..." % MinDist
    scaffolds = list(set([i[1] for i in exons_dict.keys()]))
    for scaff in scaffolds:
        genes_in_scaff = []
        for key in exons_dict:
            if key[1] == scaff:
                gene_loci = (key[0],) + tuple(sorted(exons_dict[key][8:10]))
                genes_in_scaff.append(gene_loci)
        for gene in filterByMinDist(genes_in_scaff, MinDist):
            del(exons_dict[(gene,scaff)])
   
    print "There are %s exons" % len(exons_dict);
    return exons_dict




def storeExonsInFrame(exons_dict, queries_db, out_file):
    """
    Strips the exon's ends, so that it is in frame, and then stores the results in a file.

    Strip the exon's first and last residues so that they correspond to the start
    and end of a codon, respectively. Then stores the stripped exons in a file.
    """
    
    queries_dict = SeqIO.index(queries_db, "fasta")

    print "Storing exons ..."
    exons_in_frame = []
    for exon in exons_dict.values():
        if exon[6] % 3 == 1:
            start = exon[6]-1
            end = exon[7] - exon[7]%3
            seq = queries_dict[exon[0]].seq[start:end]
            exons_in_frame.append(
                SeqRecord(seq, id=queries_dict[exon[0]].id))

        else:
            start = exon[6] + (3 - exon[6]%3)%3
            end = exon[7] - exon[7]%3
            seq = queries_dict[exon[0]].seq[start:end]
            exons_in_frame.append(
                SeqRecord(seq, id=queries_dict[exon[0]].id))

    SeqIO.write(exons_in_frame, open(out_file, "w"), "fasta");
    print "A total of %s exons are kept" % len(exons_in_frame);
    print "These exons have been stored in the file: " + out_file;



##### USE #####
# Esto debe ir en la guia rapida!.

blast_table = "pulled_seqs_blastn_out.csv"
queries_db = "pulled_seqs.fa"
out_file = "exons_in_frame.fasta"

exons_dict = getLargestExon(blast_table) # Get exon with length > 300
exons_dict = eraseFalsePosi(exons_dict) # Erase presumable false positives.
exons_dict = wellSeparatedExons(exons_dict) # Keep exons separated by > 810KB
storeExonsInFrame(exons_dict, queries_db, out_file) # Store exons in frame.
##################


### PSEUDO-CODE ###
# PARSE blast_table
#   FOR query-sbj pair:
#       keep exon whose length >= 300 and
#           <= largest_exon from query-sbj pair
#
#   FOR remaining query-sbj pairs:
#       IF a query aligns with more than 2 distinct sbj:
#           KEEP as true positive the sbj containing the
#               greatest number of aligned exons.
#
#   FOR remaining query-sbj pairs:
#       FOR all queries matched with a given sbj:
#           KEEP query if is separated from its
#               precendent query-neighbour by more than 810KB
#
#   FOR remaining query-sbj pairs:
#       CUT borders of exon so that it is in frame
#
# STORE exons' sequences in a fasta file. 
