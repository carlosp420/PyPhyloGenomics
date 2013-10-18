'''
======
MUSCLE
======

Reads in sequences from files, group homologous sequences based on gene IDs and aligns them.
'''

from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import MuscleCommandline
import subprocess
import sys
import time
import os
import string
import re
import requests
import glob


def batchAlignment(files):
    '''
    Reads in sequences from files, group homologous sequences and aligns them.

    ``files`` a list of file names. The first file in the list should
    be the master file (e.g. exon sequences of *Bombyx mori*).

    Example:
    
    >>> from pyphylogenomics import MUSCLE
    >>> files = ['Bmori_exons.fasta', 'Danaus_exons.fasta','Heliconius_exons.fasta','Manduca_exons.fasta']
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

            path = os.path.join(folder, ID1.split(':')[0] + ".fasta");
            muscle_cline = MuscleCommandline(out=path, fasta=True)
        
            child = subprocess.Popen(str(muscle_cline),
                                     stdin = subprocess.PIPE,
                                     stderr = subprocess.PIPE,
                                     shell = (sys.platform!="win32")) # In case the user's machine is not Windows.

            SeqIO.write(seqs_list, child.stdin, "fasta")
            child.stdin.close()

            time.sleep(3) # Wait 3 seconds to not fill off the screen with CMDs.

    print "%d alignments have been saved in the folder \"%s\"" % (count, folder)




def bluntSplicer(folder_path, window=20): 
    '''
    Splices Multiple Sequence Alignments objects found in *folder_path*.
    Objects should be in FASTA format. Gaps from both flanks of the alignments
    will be trimmed.

    Window size is used to find gaps in flanks (default 20 nucleotides).
    **bluntSplicer** reads the MSA files from *folder_path*, and stores the spliced MSAs
    in the same folder.

    Example:

    >>> from pyphylogenomics import MUSCLE
    >>> MUSCLE.bluntSplicer("alignments/") # folder_path containing the FASTA file alignments
    '''
    
    folder_path = folder_path.strip();
    folder_path = re.sub("/$", "", folder_path);

    # MUSCLE.batchAlignment.py created a alignments folder_path.
    for path in glob.glob(os.path.join(folder_path, "*.fasta")): 
        alignment = AlignIO.read(path,"fasta")
        print "\nSplicing %s file" % path.split("\\")[-1]
        for i in range(0,alignment.get_alignment_length()): # For checking gaps on the left flank.
            gap = False
            for record in alignment[:,i:i+window]:
                if "-" in record.seq:
                    gap = True
                    break
                else:
                    continue
            if gap == False:
                start = i
                break

        for i in range(alignment.get_alignment_length(),0,-1): # For checking gaps on the right flank.
            gap = False
            for record in alignment[:,i-window:i]:
                if "-" in record.seq:
                    gap = True
                    break
                else:
                    continue
            if gap == False:
                end = i
                break

        out = path.split(".fasta")[0] + "_bluntlySpliced.fasta" # Saving edited MSAs in same folder.
        AlignIO.write (alignment[:,start:end], out, "fasta")




def designPrimers(folder, tm="55", min_amplength="100", max_amplength="500", gencode="universal", mode="primers", clustype="dna", amptype="dna_GTRG", email=""):
    '''
    It will send a FASTA alignment to `primers4clades`_ in order to design degenerate primers. Input data needed:

    * Alignment in FASTA format containing at least 4 sequences.
    * Several parameters: 

        * temperature
        * minimium amplicon length
        * maximum amplicon length
        * genetic code
        * cluster type
        * substitution model
        * email address
        
   Example:
   The values shown are the default. Change them if needed.

    >>> from pyphylogenomics import MUSCLE

    >>> folder = "alignments"   # folder containing the FASTA file alignments
    >>> tm = "55"               # annealing temperature
    >>> min_amplength = "250"   # minimium amplicon length
    >>> max_amplength = "500"   # maximum amplicon length
    >>> gencode = "universal"   # see below for all available genetic codes
    >>> mode  = "primers"
    >>> clustype = "dna"
    >>> amptype = "dna_GTRG"    # substitution model used to estimate phylogenetic information
    >>> email = "youremail@email.com"   # primer4clades will send you an email with very detailed results

    >>> MUSCLE.designPrimers(folder, tm, min_amplength, max_amplength, gencode, mode, clustype, amptype, email)

    The best primer pairs will be printed to your screen. Detailed results will
    be saved as HTML files in your alignments folder. But it is recommended if
    you also get the results by email. primers4clades_ will send you one email
    for each alignment.

    The genetic code table (variable ``gencode``) can be any of the following:

    * ``universal`` for standard
    * ``2`` for vertebrate mitochondrial
    * ``3`` for yeast mitochondrial
    * ``4`` for mold and protozoa mitochondrial
    * ``5`` for invertebrate mitochondrial
    * ``6`` for ciliate
    * ``9`` for echinoderm and flatworm
    * ``10`` for  euplotid nuclear
    * ``11`` for  bacterial and plastid
    * ``12`` for  alternative yeast nuclear
    * ``13`` for  ascidian mitochondrial
    * ``14`` for  flatworm mitochondrial
    * ``15`` for  Blepharisma nuclear
    * ``16`` for  Chlorophycean mitochondrial
    * ``21`` for  Trematode mitochondrial
    * ``22`` for  Scenedesmus obliquus mitochondrial
    * ``23`` for  Thraustochytrium mitochondrial

    The evolutionary substitution model can be any of the following (variable
    ``amptype``):

    * ``protein_WAGG``  for protein WAG+G
    * ``protein_JTTG``  for protein JTT+G
    * ``protein_Blosum62G``  for protein Blosum62+G
    * ``protein_VTG``  for protein VT+G
    * ``protein_DayhoffG``  for protein Dayhoff+G
    * ``protein_MtREVG``  for protein MtREV+G
    * ``dna_HKYG``  for dna HKY+G
    * ``dna_GTRG``  for dna GTR+G
    * ``dna_K80G``  for dna K80+G
    * ``dna_TrNG``  for dna TrN+G
    * ``dna_JC69G``  for dna JC69+G

    .. _primers4clades: http://floresta.eead.csic.es/primers4clades/#0
    '''

    if os.path.exists(folder):
        # are there alignments/files in that folder?
        all_files = os.path.join(folder, "*")
        alns = glob.glob(all_files)

        if len(alns) > 0:
            url = "http://floresta.eead.csic.es/primers4clades/primers4clades.cgi"
            params = { 'tm': tm,  'min_amplength': min_amplength,
                    'max_amplength': max_amplength, 'mode': mode, 'gencode':
                    gencode, 'clustype': clustype, 'email': email}

            primers = [];
            for aln in alns:
                # match only files ending in .FAS[TA]
                if re.search("fas[ta]*$", aln, re.I):
                    print "\nProcessing file \"%s\"" % aln
                    files = {'sequencefile': open(aln, 'rb')}
                    r = requests.post(url, files=files, data=params);
    
                    this_file = os.path.split(aln)[1];
                    this_file = re.sub(".fas.*", "", this_file);

                    # Save result to file
                    to_print  = "Writing detailed results as file \"";
                    to_print += str(aln) + ".html\"";
                    print to_print;

                    f = open(str(aln) + ".html", "w");
                    f.write(r.text);
                    f.close()

                    # Show primer pair to user
                    html_file = string.split(r.text, "\n");
                    i = 1;
                    for line in html_file:
                        if "degen_corr" in line:
                            seq = "";
                            seq = line.split(" ")[0].strip();

                            description = line.split(" ")[2].strip();

                            this_id  = this_file + "_" + line.split(" ")[1].strip();
                            this_id += "_" + str(i);

                            seq = Seq(seq, generic_dna);
                            seq_record = SeqRecord(seq);
                            seq_record.id = this_id;
                            seq_record.description = description;
                            primers.append(seq_record);
                            i = int(i);
                            i = i + 1;

                        if i == 3:
                            break;


            # Write primers to alignment file
            SeqIO.write(primers, "primers.fasta", "fasta");
            print "\nDone.\nAll primers have been saved in the file \"primers.fasta\"";

        else:
            print "\nError! the folder \"%s\" is empty.\n" % folder;

    else:
        print "\nError! the folder \"%s\" does not exists.\n" % folder;

