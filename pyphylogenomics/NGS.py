#!/usr/bin/env python

'''
===
NGS
===

Prepares output data from sequencing round in IonTorrent (FASTQ format file). 

* Changes quality format from Phred to Solexa (which is required by the fastx-toolkit). 
* Changes sequences id to incremental numbers.
* Creates temporal FASTA file.

This module also finds matches between reads and expected gene regions by using
BLASTn, and separates reads based on matches against indexes.

It also has functions to do assembly of reads using ``velvet``.
'''

import glob;
import sys;
import os;
import re;
import shutil;
from Bio import SeqIO;
import subprocess;
import requests;
import pp;


def prepare_data(ionfile, index_length):
    '''
    * Changes quality format from Phred to Solexa (which is required by the fastx-toolkit). 
    * Changes sequences id to incremental numbers.
    * Creates temporal FASTA file with the indexes removed from the sequences.

    Files generated will be written to folder ``data/modified/`` 

    * ``ionfile`` argument is FASTQ format file as produced by IonTorrent
    * ``index_length`` number of base pairs of your indexes. This is necessary \
                       to trim the indexes before blasting the FASTA file      \
                       against the reference gene sequences.

    Example:

    >>> from pyphylogenomics import NGS
    >>> ionfile = "ionrun.fastq";
    >>> index_length = 8;
    >>> NGS.prepare_data(ionfile, index_length);
    Your file has been saved using Solexa quality format as data/modified/wrk_ionfile.fastq
    Your sequence IDs have been changed to numbers.
    The FASTA format file data/modified/wrk_ionfile.fasta has been created.
    '''
    # create folder to keep data
    folder = os.path.join("data", "modified");
    if not os.path.exists(folder):
        os.makedirs(folder);

    # change quality format from Phred to Solexa (required by fastx-toolkit)    
    # write file to work on
    wrkfile = os.path.join(folder, "wrk_ionfile.fastq")
    SeqIO.convert(ionfile, "fastq", wrkfile, "fastq-solexa");
    print "Your file has been saved using Solexa quality format as " + wrkfile

    # change sequences id to incremental numbers
    command = "fastx_renamer -n COUNT -i " + wrkfile + " -o tmp.fastq"
    p = subprocess.check_call(command, shell=True);
    if p != 0:
        print "\nError, couldn't execute " + command;
        sys.exit();
    print "Your sequence IDs have been changed to numbers."

    # replace working file with temporal file
    os.rename("tmp.fastq", wrkfile);

    # create temporal FASTA file
    command = "fastq_to_fasta -i " + wrkfile + " -o tmp.fasta";
    p = subprocess.check_call(command, shell=True);

    # trim index region
    index_length = int(index_length) + 1;
    command  = "fastx_trimmer -f " + str(index_length) + " -i tmp.fasta " 
    command += "-o " + os.path.join(folder, "wrk_ionfile.fasta");
    p = subprocess.check_call(command, shell=True);

    if os.path.isfile("tmp.fasta"):
        os.remove("tmp.fasta");

    print "The FASTA format file " + os.path.join(folder, "wrk_ionfile.fasta") \
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

    This step will accept matching reads that align more than 40bp to the
    expected gene sequence. Function :py:func:`NGS.filter_reads`

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

    # ppservers list to auto-discovery
    ppservers = ("*",)
    job_server = pp.Server(ppservers=ppservers, secret="123")
    jobs = []

    # split fastq file into chunks
    for blast_chunk in glob.iglob(os.path.join("output", "_re*.csv")):
        jobs.append(job_server.submit(split_ionfile_by_results, (ion_file, blast_chunk),(),()))

    for job in jobs:
        job();

    job_server.print_stats();

    # progressbar
    progressbar_width = 20;
    sys.stdout.write("Progress: [%s]" % (" " * progressbar_width))
    sys.stdout.flush();
    sys.stdout.write("\b" * (progressbar_width + 1))

    # now we need to iterate over each chunk of blast result and fastq file
    # and separate the reads into bins according to matching gene
    # filter_reads(ion_chunk, blast_chunk):
    # folder is the "output" folder that we are using to keep our result data
    jobs = []
    for ion_chunk in glob.iglob(os.path.join("output", "_re*.fastq")):
        blast_chunk = ion_chunk[:-5] + "csv"
        jobs.append(job_server.submit(filter_reads, (ion_chunk, blast_chunk, folder),(),()))

    for job in jobs:
        sys.stdout.write("#")
        sys.stdout.flush();
        job();

    sys.stdout.write("\n")
    job_server.print_stats();



def separate_by_index(fastq_file, index_list, folder="", levenshtein_distance=1):
    from Bio import SeqIO
    from pyphylogenomics import NGS

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
            found_index = NGS.find_index_in_seq(seq_record, fastq_record, levenshtein_distance);
            if found_index == "TRUE":
                basename = os.path.basename(fastq_file);
                if folder != "":
                    filename = "index_" + str(seq_record.id) + "_" + re.sub(".fastq", "", basename) + ".fastq";
                    filename = os.path.join(folder, filename);
                else:
                    filename = "index_" + str(seq_record.id) + "_" + re.sub(".fastq", "", basename) + ".fastq";
                output_handle = open(filename, "a");
                SeqIO.write(fastq_record, output_handle, "fastq");
                output_handle.close();


def assembly(fastq_file, index_length, min_quality=20, percentage=70, min_length=50):
    '''
    Do *de novo* assembly of expected sequences after doing quality control of a
    FASTQ file. Quality control includes dropping reads with low quality values,
    trimming of bad quality end and trimming index region.

    * ``fastq_file`` FASTQ format file that ideally has been separated by gene and index using the functions :py:func:`NGS.parse_blast_results` and :py:func:`NGS.separate_by_index`
    * ``index_length`` number of base pairs of the indexes. They will be trimmed from the reads during processing.
    * ``min_quality`` minimum quality score to keep (20 by default)
    * ``percentage``  minimum percent of base pairs that need to have ``min_quality`` (70% by default)
    * ``min_length`` minimum length of read sequences to keep (50 by default)
    

    Example:

    >>> from pyphylogenomics import NGS;
    >>> fastq_file   = "index_Ion_4_gene_rps5.fastq";
    >>> index_length = 8;
    >>> min_quality  = 30; # optional
    >>> percentage   = 80; # optional
    >>> min_length   = 60; # optional
    >>> NGS.assembly(fastq_file, index_length, min_quality, percentage, min_length);
    '''
    # Do quality control first
    q = quality_control(fastq_file, index_length, min_quality, percentage, min_length);
    try:
        with open('assembly_velvet.sh'): pass
    except:
        # downloading scripts to local folder
        r = requests.get("https://raw.github.com/carlosp420/PyPhyloGenomics/master/assembly_velvet.sh")
        f = open("assembly_velvet.sh", "w");
        f.write(r.content);
        f.close();

        r = requests.get("https://raw.github.com/carlosp420/PyPhyloGenomics/master/assembly_velvet2.sh")
        f = open("assembly_velvet2.sh", "w");
        f.write(r.content);
        f.close();

    if q == "ok":
        command = "bash assembly_velvet.sh filter3.fastq";
        filter3_output = subprocess.check_output(command, shell=True);
        filter3_params = get_velvet_params(filter3_output);

        best_input_kmer = guess_best_kmer(filter3_params);
        command = "bash assembly_velvet2.sh " + best_input_kmer[0] + ".fastq "
        command += str(best_input_kmer[1]);
        assembly = subprocess.check_call(command, shell=True);

        if assembly == 0:
            if count_reads("test/contigs.fa", "fasta") > 0:
                print "The assembly produced " + str(count_reads("test/contigs.fa", "fasta")) + " potential contigs";
                filename = re.sub(".fastq$", "", fastq_file) + "_assembled.fasta";
                os.rename("test/contigs.fa", filename);
                print "Assembled sequence has been saved as file " + filename;
    else:
        print "Couldn't process file " + fastq_file;


def count_reads(fastqFile, file_format):
    '''
    \* *Internal function* \*
    '''
    count = 0;
    for seq_record in SeqIO.parse(fastqFile, file_format):
        count = count + 1;
    return count;

# ----------------------------------------------------------------------------
# @input: output from runing velvet assembly on all Kmer values
# @output: a dictionary with the parameters: kmer, nodes, n50, max, total
def get_velvet_params(output):
    '''
    \* *Internal function* \*
    '''
    output = output.split("\n");
    mydict = dict();
    kmer = 31;
    for line in output:
        if "n50" in line:
            lista = dict();
            nodes = re.search("(\d+)\snodes", line);
            nodes = nodes.groups()[0];
            lista['nodes'] = nodes;

            n50 = re.search("n50 of (\d+)", line);
            n50 = n50.groups()[0]; 
            lista['n50'] = n50;

            maxim = re.search("max\s(\d+)", line);
            maxim = maxim.groups()[0];
            lista['max'] = maxim;

            total = re.search("total\s(\d+)", line);
            total = total.groups()[0];
            lista['total'] = total;

            mydict[kmer] = lista;
            kmer = kmer - 2;
    return mydict;
                


# ----------------------------------------------------------------------------
# @input: params from two runs of velvet on all Kmer values
#            these inputs are dictionaries
# @output: a list containing:
#            - the filtered file number either filter2 or filter3
#            - the best kmer value found by comparison of the two
def guess_best_kmer(filter3_params):
    '''
    \* *Internal function* \*
    '''
    n50 = [];
    for i in filter3_params:
        n50.append(int(filter3_params[i]['n50']));
        n50.sort();
        n50.reverse();
    filter3_n50 = n50[0];

    for i in filter3_params:
        if filter3_params[i]['n50'] == str(filter3_n50):
            return ["filter3", i];
        
# ----------------------------------------------------------------------------
# @input: output from runing velvet assembly on all Kmer values
# @output: a dictionary with the parameters: kmer, nodes, n50, max, total
def get_velvet_params(output):
    '''
    \* *Internal function* \*
    '''
    output = output.split("\n");
    mydict = dict();
    kmer = 31;
    for line in output:
        if "n50" in line:
            lista = dict();
            nodes = re.search("(\d+)\snodes", line);
            nodes = nodes.groups()[0];
            lista['nodes'] = nodes;

            n50 = re.search("n50 of (\d+)", line);
            n50 = n50.groups()[0]; 
            lista['n50'] = n50;

            maxim = re.search("max\s(\d+)", line);
            maxim = maxim.groups()[0];
            lista['max'] = maxim;

            total = re.search("total\s(\d+)", line);
            total = total.groups()[0];
            lista['total'] = total;

            mydict[kmer] = lista;
            kmer = kmer - 2;
    return mydict;
                




def quality_control(fastq_file, index_length, min_quality=20, percentage=70, min_length=50):
    '''
    \* *Internal function* \*

    Use fastx-tools to do quality control on FASTQ file.
    '''

    command  = "fastq_quality_filter -q " + str(min_quality) + " -p " + str(percentage);
    command += " -i " + fastq_file + " -o filter1.fastq";
    p = subprocess.check_call(command, shell=True);
    print "Processing file " +  fastq_file;

    index_length = index_length + 1;
    command = "fastx_trimmer -f " + str(index_length) + " -i filter1.fastq -o filter2.fastq";
    p = subprocess.check_call(command, shell=True);
    print "Removing indexes"

    command = "fastq_quality_trimmer -t " + str(min_quality) + " -l ";
    command += str(min_length) + " -i filter2.fastq -o filter3.fastq";
    p = subprocess.check_call(command, shell=True);
    print "Trimming low quality end";

    if p == 0:
        return "ok";
	


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
        barcode_read = str(seq.seq)[0:len(barcode_seq)];

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
    # _reaa.fastq
    ion_chunk = open(blast_chunk[:-4] + ".fastq", "w");

    fp = open(ion_file);
    for i, line in enumerate(fp):
        if i >= first_line:
            ion_chunk.write(line);
        if i > last_line - 1:
            print "Splitting " + ion_file + " into chunk " + blast_chunk[:-4] + ".fastq"
            break
    fp.close();

    ion_chunk.close();

    # return ion_chunk filename
    return blast_chunk + ".fastq";


    
def prune(folder, blast_data, seq_record, ion_id, min_aln_length):
    '''
    \* *Internal function* \*

    Takes a list of BLAST results and gets the gene_ids to save the current
    FASTQ seq_record into a file.
    It alse removes from the list those results that have been saved to a file
    and returns the list.
    '''
    from Bio import SeqIO;
    list = []
    for i in blast_data:
        if str(i.split(",")[0]) == str(ion_id):
            dic = {}
            dic['line'] = i
            i = i.split(",")
            dic['gene_id'] = i[1]
            dic['ion_id'] = ion_id
            dic['al_length'] = i[3]
            list.append(dic)

    # avoid saving the same seq_record into the gene bin twice 
    last_saved = ""
    for i in list:
        if last_saved != i['gene_id']:
            if i['al_length'] > min_aln_length:
                filename = os.path.join(folder, "gene_" + i['gene_id'] + ".fastq");
                exon_out = open(filename, "a");
                SeqIO.write(seq_record, exon_out, "fastq");
                exon_out.close();
                last_saved = i['gene_id']
        blast_data.remove(i['line'])

    return blast_data




def filter_reads(ion_chunk, blast_chunk, folder):
    from Bio import SeqIO;
    '''
    \* *Internal function* \*

    Accepting alignment lengths higher than 40 bp
    longer than our primer lengths
    '''
    min_aln_length = 40;

    blast_file = open(blast_chunk, "r");
    tmp = blast_file.readlines();
    blast_file.close();

    blast_data = []
    for i in tmp:
        blast_data.append(i.strip())
        

    # iterate over ion torrent reads
    for seq_record in SeqIO.parse(ion_chunk, "fastq"):
        if len(blast_data) > 0:
            #print "\n\nNew record--------------------"
            #print "seq record id @%s" % seq_record.id
            # avoid processing seq_records that are not in blast file
            # first id in blast_data
            #print blast_data
            first_id_in_blast_data = blast_data[0].split(",")[0]
            #print "fist id in blast_data %s" % first_id_in_blast_data

            if int(seq_record.id) >= int(first_id_in_blast_data):
                #if str(seq_record.id) == ion_id and aln_length > min_aln_length:
                if str(seq_record.id) == first_id_in_blast_data:
                    #print "prune"
                    blast_data = prune(folder, blast_data, seq_record,
                                first_id_in_blast_data, min_aln_length)
                else:
                    break

