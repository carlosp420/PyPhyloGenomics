#!/usr/bin/env python

'''
===
NGS
===

Prepares output data from sequencing round in IonTorrent (FASTQ format file). 
(i) Changes quality format from Phred to Solexa (which is required by the 
    fastx-toolkit). 
(ii) Changes sequences id to incremental numbers.
(iii) Creates temporal FASTA file.
'''

import glob;
import sys;
import os;
import shutil;
from Bio import SeqIO;
import subprocess;


def prepare_data(ionfile):
    '''
    * Changes quality format from Phred to Solexa (which is required by the 
      fastx-toolkit). 
    * Changes sequences id to incremental numbers.
    * Creates temporal FASTA file.
    Files generated will be written to folder ``data/modified/`` 

    ``ionfile`` argument is FASTQ format file as produced by IonTorrent

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


def parse_blast_results(blast_table, sbj_db):
    '''
    This function uses the BLAST results from a CSV file and separates the 
    IonTorrent reads into bins that match each of the target genes (or reference
    sequences).
    To speed up things a few tricks are made:
        * Split CSV file into several chunks.
        * Split ionfile.fastq into complementary chunks so that the reads in the CSV chunks
          will be found in our fastq chunks. 
    '''
    # create folder to keep output data
    folder = "output";
    if not os.path.exists(folder):
        os.makedirs(folder);

    p = split_results_file(blast_table);
    if p != 0:
        print "\nError, couldn't execute parse_blast_results";
        sys.exit();

    # move splitted blast results to output/
    for file in glob.iglob("_re*"):
        os.rename(file, os.path.join("output", file));



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
        head = head.split("\t");
        first_id = "@" + head[0];
    
    # get last id from result chunk
    number_of_lines_in_chunk = sum(1 for line in open(blast_chunk))
    fp = open(blast_chunk)
    for i, line in enumerate(fp):
        if i == number_of_lines_in_chunk - 1:
            last_id = "@" + line.split("\t")[0]
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
    fp = open(ion_file);
    for i, line in enumerate(fp):
        if i >= first_line:
            print line.strip();
        if i > last_line - 1:
            break
    fp.close();
    
