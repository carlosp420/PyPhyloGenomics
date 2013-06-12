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

    >>> from pyphylogenomics import NGS
    >>> ionfile = "ionrun.fastq";
    >>> NGS.prepare_data(ionfile);
    Your file has been saved using Solexa quality format as data/modified/wrk_ionfile.fastq
    Your sequence IDs have been changed to numbers.
    The FASTA format file data/modified/wrk_ionfile.fasta has been created.

Find reads matching target genes
--------------------------------
We can separate the sequenced reads that match the expected genes by using BLAST. For this,
we need as input a FASTA format file to create a BLAST database.
We will blast the file ``wrk_ionfile.fasta`` and then will parse the results to divide our
IonTorrent data in several bins (one bin per gene).

    >>> from pyphylogenomics import BLAST;
    >>> query_seqs = "data/modified/wrk_ionfile.fasta";
    >>> genome = "target_genes.fasta";
    >>> BLAST.blastn(query_seqs, genome); 
    BLASTn finished!
    The BLAST results were written in to the file  data/modified/wrk_ionfile_blastn_out.csv
    
    
