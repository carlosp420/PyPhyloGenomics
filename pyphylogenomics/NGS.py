#!/usr/bin/env python

'''
===
NGS
===

Prepares output data from sequencing round in IonTorrent (FASTQ format file). 
(i) Changes quality format from Phred to Solexa (which is required by the fastx-toolkit). 
(ii) Changes sequences id to incremental numbers.
(iii) Creates temporal FASTA file.
'''

import os;
import shutil;
from Bio import SeqIO;
import subprocess;


def prepare_data(ionfile):
    '''
    * Changes quality format from Phred to Solexa (which is required by the fastx-toolkit). 
    * Changes sequences id to incremental numbers.
    * Creates temporal FASTA file.
    Files generated will be written to folder ``data/modified/`` 

    ``ionfile`` argument is FASTQ format file as produced by IonTorrent
    '''
    # create folder to keep data
    folder = os.path.join("data", "modified");
    if not os.path.exists(folder):
        os.makedirs(folder);

    # change quality format from Phred to Solexa (required by fastx-toolkit    
    # write file to work on
    wrkfile = os.path.join(folder, "wrk_ionfile.fastq")
    SeqIO.convert(ionfile, "fastq", wrkfile, "fastq-solexa");

    # create temporal FASTA file
    command = "fastq_to_fasta -i " + wrkfile + " -o " + os.path.join(folder,
                "wrk_ionfile.fasta");
    p = subprocess.check_output(command, shell=True);
    print p;
