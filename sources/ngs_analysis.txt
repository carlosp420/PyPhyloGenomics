=============================================
 Analysis of raw data output from IonTorrent
=============================================

We assume that you have ordered primers and perfomed PCR reactions of the found genes on your voucher specimens. 
We also assume that you have followed any wet-lab protocol for library preparation for IonTorrent sequencing. 
In the paper describing PyPhyloGenomics (doi:XXXXXXXXXX), we have used a modified version of 
the protocol in Meyer & Kircher (2010) [1]_.

.. [1] Meyer M., Kircher M. 2010. Illumina Sequencing Library Preparation for Highly Multiplexed Target Capture and Sequencing. Cold Spring Harbor Protocols. 2010:5448.

We assume that your sequencing round was successful and you have data to analyze. In this tutorial,
we will show how we analyzed our IonTorrent sequence data using PyPhyloGenomics.

Prepare raw NGS data
--------------------
The function ``prepare_data()`` in the module ``NGS`` will make a copy of your NGS data (which should be in
FASTQ format) into a new file and do the following:

* Change the quality format from Phred to Solexa (which is required by the fastx-toolkit).
* Change the sequences id to incremental numbers.
* Create a FASTA file for temporal use.

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

Find reads matching target genes
--------------------------------
We can separate the sequenced reads that match the expected genes by using BLAST. For this,
we need as input a FASTA format file to create a BLAST database.
We will blast the file ``wrk_ionfile.fasta`` and then will parse the results to divide our
IonTorrent data in several bins (one bin per gene).

This step will accept matching reads that align more than 40bp to the
expected gene sequence. Function :py:func:`NGS.filter_reads`

    >>> from pyphylogenomics import BLAST;
    >>> query_seqs = "data/modified/wrk_ionfile.fasta";
    >>> genome = "target_genes.fasta";
    >>> BLAST.blastn(query_seqs, genome); 
    BLASTn finished!
    The BLAST results were written in to the file  data/modified/wrk_ionfile_blastn_out.csv

We will use the recently generated ``wrk_ionfile_blastn_out.csv`` file to filter our reads
that match our expected gene sequences.

    >>> from pyphylogenomics import NGS;
    >>> blast_table = "data/modified/wrk_ionfile_blastn_out";
    >>> ion_file    = "data/modified/wrk_ionfile.fastq";
    >>> NGS.parse_blast_results(blast_table, ion_file);

It will take a while to parse the results. The output will be several FASTQ files (one
per target gene) containing our matching IonTorrent reads. The files will have the 
prefix **gene_**.


Separate gene bins according to indexes (or barcodes)
-----------------------------------------------------
Then we need to separate the files prefixed by **gene_** according to the indexes
that were used in our wet-lab protocol.
We used one index (or barcode) for each voucher specimen that went into the 
Ion Torrent.
And we will use those reads to figure out which reads belong to each of our specimens.
It is possible that due to some errors in base calling during the sequencing procedure 
mistakes might appear in the sequenced index region. 
Thus, we will need to perform string comparisons accepting differences of up to 1 nucleotide
between our expected and sequenced indexes.
Our indexes differ in two nucleotides (they should ideally differ more) so it is safe
to accept up to 1 mistake during the sequencing of the index region.

PyPhyloGenomics uses `Levenshtein distances <http://en.wikipedia.org/wiki/Levenshtein_distance>`_ 
for comparison of index sequences.

We assume that our FASTQ bins to separate are in the folder ``output`` and begin with the 
prefix ``gene``.
    
    >>> from pyphylogenomics import NGS;
    >>> import glob; # this module allow us selecting many files by using wildcards
    >>> index_list           = "indexes.fasta";
    >>> folder               = "output";
    >>> levenshtein_distance = 1; # maximum number of differences allowed between index comparisons
    >>> for file in glob.glob("output/gene*.fastq"):
    ...     NGS.separate_by_index(file, index_list, folder, levenshtein_distance);

All the reads will be saved to files with the prefix **index_** and a identifier
depending on the sequence ID found in the ``indexes.fasta`` file.

For example, one of our files was written as **index_IonADA_2_gene_rps5.fastq**.

Depending on the amount of your data, this process can take several hours.


Assembly of reads into consensus sequences
------------------------------------------
Once we have one FASTQ file for each gene and for each specimen or barcode, we can
start doing the assembly of the consensus sequence.

It might be necessary for you to try out several parameters that will affect the final
assembled sequence. For example you might want to be strict and remove read with 
low quality values (by default PyPhyloGenomics uses ``min_quality = 20``) or change
the percentage of high quality values per read (default is ``percentage = 80``).

After each iteration of parameters you can check the assembled sequences and 
evaluate the sequence length and coverage and decide whether to accept of reject 
the sequence.

We will do the assembly using the commonly used assembler ``velvet``:

    >>> from pyphylogenomics import NGS;
    >>> fastq_file   = "index_Ion_4_gene_rps5.fastq";
    >>> index_length = 8;
    >>> min_quality  = 30; # optional
    >>> percentage   = 80; # optional
    >>> min_length   = 60; # optional
    >>> NGS.assembly(fastq_file, index_length, min_quality, percentage, min_length);
    Final graph has 3 nodes and n50 of 102, max 121, total 250, using 0/60354 reads
    The assembly produced 2 potential contigs
    Assembled sequence has been saved as file index_Ion_4_gene_rps5_assembled.fasta

This procedure performs quality control over the reads and then trying out 
assembling the sequences using several Kmer values. Later it tries to guess the 
best Kmer value and use it to do a final assembly. 

The potential contigs are written to the file  ``index_Ion_4_gene_rps5_assembled.fasta``.

