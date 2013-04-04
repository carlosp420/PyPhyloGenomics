Getting started with PyPhyloGenomics
====================================

Some snippets of code to get you started with writing code using PyPhylogenomics.

We need to obtain candidate genes to be used in phylogenetic inference that have to fulfill the following requirements:

* Our genes should be orthologs.
* Our genes should be single-copy genes.
* Their sequence need to be around 251 DNA base pairs in length.

We will assume that our Next Generation Sequencer available is the IonTorrent_.

We have to consider the IonTorrent_ platform requirements to arrive to our target 250bp gene length:

====================  ===========
Primer                Length (bp)
====================  ===========
Adapter A             30
5' Index              8
5' Degenerate Primer  25
Exon                  **???**
3' Degenerate Primer  25
3' Index              8
Adapter P             23
====================  ===========

For IonTorrent_ Platform2, the maximum length that can be sequenced is from 280bp to 320bp in total. Thus, ``320 - 119 = 201`` is the maximum internal gene region (region within degenerate primers).

Therefore, for the new set of primers, being designed for Platform2, we have a maximum amplicon size of ``201 + 25*2 = 251bp``. 

The OrthoDB_ database has a catalog of orthologous protein-coding genes for vertebrates, arthropods and other living groups.

.. _IonTorrent: http://www.iontorrent.com/
.. _OrthoDB: http://cegg.unige.ch/orthodb6
.. _OrthoDB6_Arthropoda_tabtex.gz: ftp://cegg.unige.ch/OrthoDB6/

For this guick getting-started guide, we will download the table of orthologous genes for Arthropoda from OrthoDB's ftp server OrthoDB6_Arthropoda_tabtex.gz_.

Starting off:

    >>> import pyphylogenomics
    >>> from pyphylogenomics import OrthoDB

To see the full documentation you can use:

    >>> help(pyphylogenomics) # shows list of modules
    >>> help(pyphylogenomics.OrthoDB) # shows help for a specific module and its functions

This also works:

    >>> help(OrthoDB) # shows help for a specific module and its functions

We will find all single-copy genes for the silk moth *Bombyx mori* using the table from OrthoDB_ as input file:

    >>> in_file = 'OrthoDB6_Arthropoda_tabtext.csv'
    >>> genes = OrthoDB.single_copy_genes(in_file, 'Bombyx mori')
    ...
    Found gene: BGIBMGA011628
    Found gene: BGIBMGA002142
    Found gene: BGIBMGA014116
    Found 12167 genes.

The variable ``genes`` is a list of IDs for single-copy genes extracted from the OrthoDB table:

    >>> genes;
    ...
    'BGIBMGA012533', 'BGIBMGA014144', 'BGIBMGA011628', 'BGIBMGA002142', 'BGIBMGA014116',
    'BGIBMGA000220', 'BGIBMGA007580', 'BGIBMGA013398', 'BGIBMGA011517', 'BGIBMGA011820',
    'BGIBMGA014318', 'BGIBMGA013470', 'BGIBMGA011304', 'BGIBMGA005643', 'BGIBMGA002698']
    >>> len(genes);
    12167

We will use these gene IDs to obtain the gene sequences from the *Bombyx mori* genome. The genome can be downloaded from silkdb.org_.
We will download the **Consensus gene set by merging all the gene sets using GLEAN: Fasta** file ``silkworm_glean_cds.fa.tar.gz``.  

Untar the gzipped file and you will get the FASTA formatted file ``silkcds.fa`` containing gene IDs and sequences for *Bombyx mori*.

.. _silkdb.org: http://www.silkdb.org/silkdb/doc/download.html

We will need to BLASTn these single-copy genes against the *Bombyx mori* genome
in order to get exon sizes:

1. Pull all sequences for our gene IDs from the CDS file and write them to a file ``pulled_seqs.fa``:

    >>> from pyphylogenomics import BLAST
    >>> cds_file = "silkcds.fa"
    >>> BLAST.get_cds(genes, cds_file)
    12167  sequences were written to file pulled_seqs.fa

2. Download the *Bombyx mori* genome from silkdb.org_ (download the file ``silkworm_genome_v2.0.fa.tar.gz``). Unzip and untar the file to your working directory and you will get the file ``silkgenome.fa``
 
3. Do a BLASTn of the sequences against the *Bombyx mori* genome. The input arguments are your file containing the sequences for single-copy genes (``pulled_seqs.fa``) and your file with the genome of *Bombyx mori* which is in FASTA format (``silkgenome.fa``).

    >>> BLAST.blastn('pulled_seqs.fa', 'silkgenome.fa')
    ...
    BLASTn finished!
    The BLAST results were written in to the file  pulled_seqs_blastn_out.csv  




