'''
======
MUSCLE
======
Reads in sequences from files, group homologous sequences based on gene IDs and aligns them.
'''

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import MuscleCommandline
import subprocess
import sys
import time
import os
import string
import re

def batchAlignment(files):
    '''
    Reads in sequences from files, group homologous sequences and aligns them.

    ``files`` a list of file names. The first file in the list should
    be the master file (e.g. exon sequences of *Bombyx mori*).

    Example:
    
    >>> from PyPhyloGenomics import MUSCLE
    >>> files = ['Bmori_exons.fas', 'Danaus_exons.fas','Heliconius_exons.fas','Manduca_exons.fas']
    >>> MUSCLE.batchAlignment(files)

    All aligned sequences will be written into a folder called ``alignments`` as
    FASTA files (one file per exon).
    '''
    
    # create folder for saving alignments
    folder = "alignments";
    if not os.path.exists(folder):
        os.makedirs(folder);


    for idx, f in enumerate(files): # Read each file into a dictionary.
        exec("handle%d = %s" % (idx+1, "SeqIO.index(f,'fasta')")) 

    count = 0 # To count the number of homologous groups.
    print "If you are using Windows..."
    time.sleep(1);
    print "many command line shells might pop up now."
    time.sleep(3)

    for ID1 in handle1: # Parsing the master file.
        seqs_list = []
        seqs_list.append(SeqRecord(handle1[ID1].seq, id = ID1))
        print "Pooling gene %s" % ID1;

        idx = len(files)
        while idx > 1: # Parsing the other files.
            script = '''for ID2 in handle%d:
    if ID1.split(':')[0] == re.sub(":.+$", "", ID2.split('-')[0]):
        seqs_list.append(SeqRecord(handle%d[ID2].seq, id = ID2))'''
            exec(script % (idx, idx))
            idx -= 1

        #TODO: Add a command for PSAs, for cases where the homologous group
                                                    # has only 2 sequences.

        if len(seqs_list) > 2: # Now perfom the MSA.
            count += 1

            path = os.path.join(folder, ID1.split(':')[0] + ".fas");
            muscle_cline = MuscleCommandline(out=path, fasta=True)
        
            child = subprocess.Popen(str(muscle_cline),
                                     stdin = subprocess.PIPE,
                                     stderr = subprocess.PIPE,
                                     shell = (sys.platform!="win32")) # In case the user's machine is not Windows.

            SeqIO.write(seqs_list, child.stdin, "fasta")
            child.stdin.close()

            time.sleep(3) # Wait 3 seconds to not fill off the screen with CMDs.

    print "%d alignments have been saved in the folder \"%s\"" % (count, folder)


