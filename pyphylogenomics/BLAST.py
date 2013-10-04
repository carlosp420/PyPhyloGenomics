'''
=====
BLAST
=====

Prepares data and executes BLAST commands to create BLAST database.
Performs blastn of user sequences against genomic sequences.
Provides functions to parse the generated BLAST tables and extract exons.
'''

import re;
import sys;
import subprocess;
from Bio import SeqIO;
from Bio.SeqRecord import SeqRecord;
from operator import itemgetter;
import pp;
import shutil;
import os;
import time;


def get_cds(genes, cds_file):
    '''
    Writes a FASTA file containing CDS sequences: ``pulled_seqs.fasta``
    
    ``genes`` argument is a list of gene IDs.
    ``cds_file`` is the file name of predicted CDS for a species in FASTA format.

    '''

    records = []; # To store the sequences that matched.
    for seq_record in SeqIO.parse(cds_file, "fasta"):
        this_id = re.sub("-TA$", "", seq_record.id);
        if this_id in genes:
            records.append(SeqRecord(seq_record.seq, id=this_id));

    SeqIO.write(records, open("pulled_seqs.fasta", "w"), "fasta") 
    print len(records), " sequences were written to file pulled_seqs.fasta in the current working directory.";



def makeblastdb(genome, mask=False):
    '''
    Creates a BLAST database from a genome in FASTA format and
    optionally eliminates low-complexity regions from the sequences.
    '''

    if mask == True:
        command = 'dustmasker -in '+ genome + ' -infmt fasta '
        command += '-outfmt maskinfo_asn1_bin -out ' + genome + '_dust.asnb'
        try:
            with open(genome + "_dust.asnb"):
                print "Using existing ``" + genome + "_dust.asnb`` database"
        except IOError:
            print "masking low_complexity regions..."
            p = subprocess.check_output(command, shell=True) # identifying low-complexity regions.
            print p
        
            command = 'makeblastdb -in ' + genome + ' -input_type fasta -dbtype nucl '
            command += '-mask_data ' + genome + '_dust.asnb '
            command += '-out ' + genome + ' -title "Whole Genome without low-complexity regions"'
            print "creating database..."
            p = subprocess.check_output(command, shell=True) # Overwriting the genome file.
            print p

    else:
        command  = 'makeblastdb -in ' + genome + ' -input_type fasta -dbtype nucl '
        command += '-out ' + genome + ' -title "Whole Genome unmasked"'
        print "creating database..."
        p = subprocess.check_output(command, shell=True) 
        print p



def blastn(query_seqs, genome, e_value=0.00001, mask=True):
    '''
    Performs a BLASTn of user's sequences against a genome. It will create a
    BLAST database from the genome file first.

    ``query_seqs`` argument is a FASTA file of users sequences.
    ``genome`` argument is a FASTA file of a species genome.
    '''

    # do a BLAST database first
    makeblastdb(genome, mask=mask);
    blast_out = query_seqs.split(".")[0]+"_blastn_out.csv" # Name of the output file.

    output  = "blasting: " + query_seqs.split("\\")[-1] + " "
    output += "against db = " + genome.split("\\")[-1]
    print output;

    # how many sequences in input file?
    nseqs = 0;
    for seq_record in SeqIO.parse(query_seqs, "fasta"):
        nseqs +=1;

    # ppservers list to auto-discovery
    ppservers = ("*",);
    job_server = pp.Server(ppservers=ppservers, secret="123");
    jobs = [];

    # split fasta file for parallel blast
    files = [];
    divisor = 1;
    if nseqs > 399:
        while nseqs/divisor > 200:
            divisor += 1;
            if divisor > 400:
                break;

        record_iter = SeqIO.parse(open(query_seqs), "fasta");

        for i, batch in enumerate(batch_iterator(record_iter, nseqs/divisor)):
            filename = "group_%i.fasta" % (i+1);
            files.append(filename);
            handle = open(filename, "w");
            count = SeqIO.write(batch, handle, "fasta");
            handle.close();
    else:
        files.append(query_seqs)


    for f in files:
        command = ('blastn -query ' + f + ' -db ' + genome + ' -task blastn -db_soft_mask 11 ' if mask 
        else 'blastn -query ' + f + ' -db ' + genome + ' -task blastn ')
        command += '-evalue ' + str(e_value) + ' -out ' + f + "_out.csv" + ' -num_threads 1 -outfmt 10'
        jobs.append(job_server.submit(do_blast, (command,), modules=('subprocess',)))

    print "\nBlasting sequences. This might take several minutes ...";

    # progressbar
    progressbar_width = int(divisor) + 2;
    sys.stdout.write("Progress: [%s]" % (" " * progressbar_width))
    sys.stdout.flush();
    sys.stdout.write("\b" * (progressbar_width + 1))

    # execute parallel jobs
    for job in jobs:
        job();
        # update the bar
        sys.stdout.write("#")
        sys.stdout.flush()

    sys.stdout.write("\n");
    
    job_server.print_stats();
    print "\nBLASTn finished!"

    # merging blast tables generated above
    destination = open(blast_out,'wb')
    for f in files:
        shutil.copyfileobj(open(f + "_out.csv",'rb'), destination)
        os.remove(f + "_out.csv")
        if nseqs > 399:
            os.remove(f)

    destination.close()
    print "The BLAST results were written into the file ", blast_out
     


def do_blast(command):
    subprocess.check_output(command, shell=True);



def getLargestExon(blast_table, E_value = 0.001, ident = 98, exon_len = 300):
    """
    Returns the highest alignment number for each query-sbj match from a blast table. 
    """
    table = open(blast_table, "r") 
    exons = {} # A dictionary, where the key is the (query,sbj) tuple
                # and the value is the longest alignment of all the possible
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
    Keeps the query-sbj match with the longest alignment number.
    
    From a dictionary generated with the getLargestExon function, where 2 or more
    query-sbj matches shared the same query, eraseFalsePosi keeps the query-sbj
    match with the longest alignment.
    """

    print "Erasing False Positives ..."
    new_exons_dict = {}

    # First sorts alphabetically by query_name and then by the total number of
    # alignments of each query-sbj match.
    keys = sorted([j + (exons_dict[j][-1],) for j in exons_dict.keys()],
                  key=itemgetter(0,2), reverse=True) 

    # Once sorted, throws the number of alignments of each query-sbj match.
    keys = [j[:2] for j in keys] 
    unique = list(set([i[0] for i in keys]))
    for i in unique:
        # Since in 'keys', the query-sbj with the highest alignment number
        # always appear on the most-left side, only this one is kept.
        for key in keys: 
            if key[0] == i:
                new_exons_dict.update({key:exons_dict[key]})
                break # That's why we use inmediately break.
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
        if loci[1] - last_end < MinDist: # Exons whose separation is less than MinDist.
            genes.append(loci[0])
            last_end = loci[2]

    return genes # This set will be filtered out in the next function.
                                              


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
        for gene in filterByMinDist(genes_in_scaff, MinDist): # Use the function above to identify those exons which don't pass the MinDist threshold.
            del(exons_dict[(gene,scaff)]) # And filter them out!
   
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
            ID = queries_dict[exon[0]].id + ':' + str(start+1)+'-'+str(end) 
            exons_in_frame.append(
                SeqRecord(seq, id=ID))

        else:
            start = exon[6] + (3 - exon[6]%3)%3
            end = exon[7] - exon[7]%3
            seq = queries_dict[exon[0]].seq[start:end]
            ID = queries_dict[exon[0]].id + ':' + str(start+1)+'-'+str(end) 
            exons_in_frame.append(
                SeqRecord(seq, id=ID))

    SeqIO.write(exons_in_frame, open(out_file, "w"), "fasta")
    print "A total of %s exons are kept" % len(exons_in_frame)
    print "These exons have been stored in the file: " + out_file



def blastParser(blast_table, sbj_db, out_file, sp_name = 'homologous', E_value = 0.01, ident = 75, exon_len = 300):
    """
    Returns the subjects' sequences aligned with the queries as long as they pass
    the thresholds given on the parameters. 

    Example:

    >>> BLAST.blastParser("LongExons_out_blastn_out.csv", "Dp_genome_v2.fasta", "Danaus_exons.fasta", sp_name="Danaus");
    Reading files ...
    Parsing BLAST table ...
    A total of 158 sequences passed the thresholds.
    They have been stored in the file: Danaus_exons.fasta

    The parameter ``sp_name`` is important as it will be used as part of the exons IDs.
    """

    print "Reading files ..."
    
    table = open(blast_table, "r")
    sbj_dict = SeqIO.index(sbj_db, "fasta")
    seqs = {}
    
    print "Parsing BLAST table ..."  
    for alignment in table:
        if alignment.startswith('Query'): continue # in case the input table has header.  

        align_vars = alignment.strip().split(",") # "blast_table" must be a csv file.
        align_vars[3:10] = [int(i) for i in align_vars[3:10]]
        align_vars[3] -= align_vars[5] # replaces the alignment length for the exon length

        if float(align_vars[10]) <= E_value and float(align_vars[2]) >= ident and align_vars[3] >= exon_len:
            if align_vars[8] < align_vars[9]: # The alignment is with the sbj positive strand.
                start = align_vars[8] - 1
                end = align_vars[9]
                seq = sbj_dict[align_vars[1]].seq[start:end]
                ID = align_vars[0] + '-' + '_'.join(sp_name.split()) + '_' + align_vars[1] + ':' + str(start+1)+'-'+str(end)
                DE = sp_name + ' sequence homologous to ' + align_vars[0] + ':' + str(align_vars[6])+'-'+str(align_vars[7])
                seqs[ID] = SeqRecord(seq, id=ID, description=DE)
                  
            else: # The alignment is with the sbj negative strand.
                start = align_vars[9] - 1
                end = align_vars[8]
                seq = sbj_dict[align_vars[1]].seq[start:end].reverse_complement() # Reverse and complement the plus strand!
                ID = align_vars[0] + '-' + '_'.join(sp_name.split()) + '_' + align_vars[1] + ':c' + str(end)+'-'+str(start+1)
                DE = sp_name + ' sequence homologous to ' + align_vars[0] + ':' + str(align_vars[6])+'-'+str(align_vars[7])
                seqs[ID] = SeqRecord(seq, id=ID, description=DE)


    SeqIO.write(seqs.values(), open(out_file, "w"), "fasta")
    print "A total of %s sequences passed the thresholds." % len(seqs)
    print "They have been stored in the file: " + out_file


def batch_iterator(iterator, batch_size) :
    """Returns lists of length batch_size.
 
    # Taken from http://biopython.org/wiki/Split_large_file

    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.
 
    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    entry = True #Make sure we loop once
    while entry :
        batch = []
        while len(batch) < batch_size :
            try :
                entry = iterator.next()
            except StopIteration :
                entry = None
            if entry is None :
                #End of file
                break
            batch.append(entry)
        if batch :
            yield batch
