#!/usr/bin/env python

'''
===
NGS
===

Prepares output data from sequencing round in IonTorrent (FASTQ format file). 

* Changes quality format from Phred to Solexa (which is required by the fastx-toolkit). 
* Changes sequences id to incremental numbers.
* Creates temporal FASTA file.
'''

import glob;
import sys;
import os;
import re;
import shutil;
from Bio import SeqIO;
import subprocess;


def prepare_data(ionfile):
    '''
    * Changes quality format from Phred to Solexa (which is required by the fastx-toolkit). 
    * Changes sequences id to incremental numbers.
    * Creates temporal FASTA file.

    Files generated will be written to folder ``data/modified/`` 

    * ``ionfile`` argument is FASTQ format file as produced by IonTorrent

    Example:

    >>> from pyphylogenomics import NGS
    >>> ionfile = "ionrun.fastq";
    >>> NGS.prepare_data(ionfile);
    '''
    # create folder to keep data
    folder = os.path.join("data", "modified");
    if not os.path.exists(folder):
        os.makedirs(folder);

    # change quality format from Phred to Solexa (required by fastx-toolkit)    
    # write file to work on
    wrkfile = os.path.join(folder, "wrk_ionfile.fastq")
    SeqIO.convert(ionfile, "fastq", wrkfile, "fastq-solexa");
    print "\nYour file has been saved using Solexa quality format as " + wrkfile

    # change sequences id to incremental numbers
    command = "fastx_renamer -n COUNT -i " + wrkfile + " -o tmp.fastq"
    p = subprocess.check_call(command, shell=True);
    if p != 0:
        print "\nError, couldn't execute " + command;
        sys.exit();
    print "\nYour sequence IDs have been changed to numbers."

    # replace working file with temporal file
    os.rename("tmp.fastq", wrkfile);

    # create temporal FASTA file
    command = "fastq_to_fasta -i " + wrkfile + " -o " + os.path.join(folder,
                "wrk_ionfile.fasta");
    p = subprocess.check_call(command, shell=True);
    print "\nThe FASTA format file " + os.path.join(folder, "wrk_ionfile.fasta") \
            + " has been created.";


def parse_blast_results(blast_table, ion_file):
    '''
    This function uses the BLAST results from a CSV file and separates the 
    IonTorrent reads into bins that match each of the target genes (or reference
    sequences).
    To speed up things a few tricks are made:
    * Split CSV file into several chunks.
    * Split ionfile.fastq into complementary chunks so that the reads in the CSV chunks will be found in our fastq chunks. 
    Reads will be written into separate FASTQ files, one per matching target
    gene.
    These files will be written in the folder ``output``.
    '''
    # create folder to keep output data
    folder = "output";
    if not os.path.exists(folder):
        os.makedirs(folder);

    # split results file into chunks
    p = split_results_file(blast_table);
    if p != 0:
        print "\nError, couldn't execute parse_blast_results";
        sys.exit();

    # move splitted blast results to output/
    for file in glob.iglob("_re*"):
        src = file;
        file += ".csv";
        dest = os.path.join("output", file); 
        os.rename(src, dest);

    # split fastq file into chunks
    for file in glob.iglob(os.path.join("output", "_re*.csv")):
        ion_chunk = split_ionfile_by_results(ion_file, file);

        # now we need to iterate over each chunk of blast result and fastq file
        # and separate the reads into bins according to matching gene
        # filter_reads(ion_chunk, blast_chunk):
        # folder is the "output" folder that we are using to keep our result data
        filter_reads(ion_chunk, file, folder);
        print "Filtering reads in file " + ion_chunk


def separate_by_index(fastq_file, index_list, folder="", levenshtein_distance=1):
    '''
    This function divides FASTQ reads into bins according to a list of indexes
    (or barcodes).
    The *index_list* should be in FASTA format.
    It will compare the template indexes and those in the reads and accept
    indexes with a difference no bigger than the *levenshtein* distance (default
    1 base pair difference).

    See http://en.wikipedia.org/wiki/Levenshtein_distance

    * ``fastq_file`` FASTQ format containing reads as produced by IonTorrent
    * ``index_list`` FASTA format file containing indexes (or barcodes)
    * ``folder`` *Optional*: Directory containing FASTQ format files to process
    * ``levenshtein_distance`` *Optional*, default = 1: Maximum number of different nucleotides that will be accepted when comparing template and sequenced indexes (due to erros in base calling during sequencing). 

    Example:

    >>> from pyphylogenomics import NGS;
    >>> fastq_file = "gene_rps5.fastq";
    >>> index_list = "indexes.fasta";
    >>> folder     = "output";
    >>> NGS.separate_by_index(fastq_file, index_list, folder);

    You can also automate parsing many FASTQ files at once:

    >>> from pyphylogenomics import NGS;
    >>> import glob; # this module allow us selecting many files by using wildcards
    >>> index_list = "indexes.fasta";
    >>> folder     = "output";
    >>> for file in glob.glob("output/gene*.fastq"):
    ...     NGS.separate_by_index(file, index_list, folder);
    '''
    print "Processing file " + fastq_file;
    if folder != "":
        folder = re.sub("/$", "", folder);
        folder = os.path.abspath(folder);
        print "Output files will be written into " + folder

    for seq_record in SeqIO.parse(index_list, "fasta"):
        for fastq_record in SeqIO.parse(fastq_file, "fastq"):
            found_index = "";
            found_index = find_index_in_seq(seq_record, fastq_record, levenshtein_distance);
            if found_index == "TRUE":
                basename = os.path.basename(fastq_file);
                if folder != "":
                    filename = str(seq_record.id) + "_" + re.sub(".fastq", "", basename) + ".fastq";
                    filename = os.path.join(folder, filename);
                else:
                    filename = str(seq_record.id) + "_" + re.sub(".fastq", "", basename) + ".fastq";
                output_handle = open(filename, "a");
                SeqIO.write(fastq_record, output_handle, "fastq");
                output_handle.close();


def find_index_in_seq(barcode, seq, levenshtein_distance):
    '''
    \* *Internal function* \*

    Arguments:  *barcode*,  *read sequence*.
    Iterate through primer's degenerated IUPAC and try to find it in sequence
    read.
    Return TRUE on success.
    '''
    found = "false";

    while (found != "true" ):
        barcode_seq = str(barcode.seq) 
        barcode_read = str(seq.seq)[0:8];

        # compare levenshtein
        distance = levenshtein(barcode_seq.upper(), barcode_read.upper());

        if distance < levenshtein_distance + 1:
            # accept
            found = "true";

            #if str(seq.id) not in fastq_id_list:
                #fastq_id_list.append(str(seq.id));
            return "TRUE";
        else:
            # reversecomplement
            barcode_seq = str(barcode.seq.reverse_complement())

            # compare levenshtein
            distance = levenshtein(barcode_seq.upper(), barcode_read.upper());

            if distance < levenshtein_distance + 1:
                #print "Found reverse complement";
                # accept
                found = "true";
    
                #if str(seq.id) not in fastq_id_list:
                    #fastq_id_list.append(str(seq.id));
                return "TRUE";

        return "FALSE";



def levenshtein(a,b):
    '''
    \* *Internal function* \*

    Calculates the Levenshtein distance between a and b.
    '''
    n, m = len(a), len(b)
    if n > m:
        # Make sure n <= m, to use O(min(n,m)) space
        a,b = b,a
        n,m = m,n
        
    current = range(n+1)
    for i in range(1,m+1):
        previous, current = current, [i]+[0]*n
        for j in range(1,n+1):
            add, delete = previous[j]+1, current[j-1]+1
            change = previous[j-1]
            if a[j-1] != b[i-1]:
                change = change + 1
            current[j] = min(add, delete, change)
            
    return current[n]



def split_results_file(blast_table):
    '''
    \* *Internal function* \*

    Split a CSV format result file from a BLAST in chunks.
    '''
    # build an slicer for Python 
    command = "split -l 10000 " + blast_table + " _re"
    p = subprocess.check_call(command, shell=True);
    return p;


def split_ionfile_by_results(ion_file, blast_chunk):
    '''
    \* *Internal function* \*

    This function divides an IonTorrent run file into chunks that match chunks
    of a BLAST results file in CSV format.
    '''

    # get first id from result chunk
    with open(blast_chunk) as myfile:
        head = myfile.readline();
        head = head.split(",");
        first_id = "@" + head[0];
    
    # get last id from result chunk
    number_of_lines_in_chunk = sum(1 for line in open(blast_chunk))
    fp = open(blast_chunk)
    for i, line in enumerate(fp):
        if i == number_of_lines_in_chunk - 1:
            last_id = "@" + line.split(",")[0]
    fp.close();
    
    # get line with first and last ids in ionfile
    fp = open(ion_file);
    first_line = 0;
    last_line  = 0;
    for i, line in enumerate(fp):
        if line.strip() == first_id:
            first_line = i;
        if line.strip() == last_id:
            # adding 3 lines due to this being FASTQ file
            last_line = i + 3;
    fp.close();
    
    # cut ionfile using first and last lines
    ion_chunk = open(blast_chunk + ".fastq", "w");

    fp = open(ion_file);
    for i, line in enumerate(fp):
        if i >= first_line:
            ion_chunk.write(line);
        if i > last_line - 1:
            print "Splitting " + ion_file + " into chunk " + blast_chunk + ".fastq"
            break
    fp.close();

    ion_chunk.close();

    # return ion_chunk filename
    return blast_chunk + ".fastq";
    


def filter_reads(ion_chunk, blast_chunk, folder):
    # get dict of ion_id : gene_id
    # to blast_data
    blast_file = open(blast_chunk, "r");
    blast_data = blast_file.readlines();
    blast_file.close();

    # iterate over ion torrent reads
    for seq_record in SeqIO.parse(ion_chunk, "fastq"):
        for line in blast_data:
            #print len(blast_data);
            orig_line = line;
            line = line.strip();
            line = line.split(",");

            ion_id  = line[0];
            gene_id = line[1];
            #print ion_id, gene_id

            if str(seq_record.id) == ion_id:
                blast_data.remove(orig_line);
                
                filename = os.path.join(folder, "gene_" + gene_id + ".fastq");
                exon_out = open(filename, "a");
                SeqIO.write(seq_record, exon_out, "fastq");
                exon_out.close();
                break;
