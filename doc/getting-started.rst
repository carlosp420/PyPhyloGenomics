Getting started with PyPhyloGenomics-python
===========================================

Some snippets of code to get you started with writing code against PyPhylogenomics.

We need to obtain candidate genes to be used in phylogenetic inference that have to fulfill the following requirements:

* Our genes should be orthologs.
* Our genes should be single-copy genes.
* Their sequence need to be around DNA 251 base pairs in length.

We will assume that our Next Generation Sequencer available is the IonTorrent_.

We have to consider the IonTorrent_ platform requirements to arrive to our target 250bp gene length:

====================  ===========
Primer                Length (bp)
====================  ===========
Adapter A             30
5' Index              8
5' Degenerate Primer  25
Exon                  ???
3' Degenerate Primer  25
3' Index              8
Adapter P             23
====================  ===========

For IonTorrent_ Platform2, the maximum length that can be sequenced is from 280bp to 320bp in total. Thus, ``320 - 119 = 201`` is the maximum internal gene region (region within degenerate primers).

Therefore, for the new set of primers, being designed for Platform2, we have a maximum amplicon size of ``201 + 25*2 = 251bp``. 

The OrthoDB_ database has a catalog of orthologous protein-coding genes for vertebrates, arthropods and other living groups.

.. _IonTorrent: http://www.iontorrent.com/
.. _OrthoDB: http://cegg.unige.ch/orthodb6
.. _OrthoDB6_Arthropoda_tabtex.gz: ftp://cegg.unige.ch/OrthoDB6/OrthoDB6_Arthropoda_tabtext.gz

For this guick getting-started guide, we will download the table of orthologous genes for Arthropoda from OrthoDB's ftp server OrthoDB6_Arthropoda_tabtex.gz_.

Starting off:

    >>> import pyphylogenomics as phylo

To see the full documentation you can use:

    >>> help(phylo) # shows list of modules
    >>> help(phylo.OrthoDB) # shows help for a specific module and its functions

We will find all single-copy genes for the silk moth *Bombyx mori* using the table from OrthoDB_ as input file.

    >>> in_file = 'OrthoDB6_Arthropoda_tabtext.csv'
    >>> genes = phylo.OrthoDB.single_copy_genes(in_file, 'Bombyx mori')



