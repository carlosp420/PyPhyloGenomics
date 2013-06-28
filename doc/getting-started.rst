===================================================
 Getting started with PyPhyloGenomics: Gene Search
===================================================

Some snippets of code to get you started with writing code using PyPhyloGenomics.

------------------------------------------
Finding candidate genes from *Bombyx mori*
------------------------------------------

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

For this quick getting-started guide, we will download the table of orthologous genes for Arthropoda from OrthoDB's ftp server OrthoDB6_Arthropoda_tabtex.gz_.

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

1. Pull all sequences for our gene IDs from the CDS file and write them to a file ``pulled_seqs.fasta``:

    >>> from pyphylogenomics import BLAST
    >>> cds_file = "silkcds.fa"
    >>> BLAST.get_cds(genes, cds_file)
    12167  sequences were written to file pulled_seqs.fasta

2. Download the *Bombyx mori* genome from silkdb.org_ (download the file ``silkworm_genome_v2.0.fa.tar.gz``). Unzip and untar the file to your working directory and you will get the file ``silkgenome.fa``
 
3. Do a BLASTn of the sequences against the *Bombyx mori* genome. The input arguments are your file containing the sequences for single-copy genes (``pulled_seqs.fasta``) and your file with the genome of *Bombyx mori* which is in FASTA format (``silkgenome.fa``).

    >>> BLAST.blastn('pulled_seqs.fasta', 'silkgenome.fa')
    ...
    BLASTn finished!
    The BLAST results were written in to the file  pulled_seqs_blastn_out.csv  

The file ``pulled_seqs_blastn_out.csv`` contains a BLAST output table with the blast results. **PyPhyloGenomics** has functions to filter out the table and get information to answer the following:

    * Which are the longest gene-sequence to genome matches?
    * Which are the genes that are "distantly enough" from each other? So that, they can used as independent evolutionary entities?

4. As stated before, we prefer long exons for each of the candidate genes ( > 300 nucleotides):

    >>> exons = BLAST.getLargestExon("pulled_seqs_blastn_out.csv", E_value=0.001, ident=98, exon_len=300)
    Parsing BLAST table ...
    Deleting exons below 300 nucleotides ...
    There are 7411 exons

5. Some small segments of sequences might match exons present in more than one scaffold. We will use the function ``eraseFalsePosi`` to keep those matches of longest length:

    >>> exons = BLAST.eraseFalsePosi(exons) # Drop presumable false positives.
    Erasing False Positives ...
    There are 6346 exons

6. Ideally we want exons that are not too close to each other in the genome to avoid gene linkage. So we will keep only those exons that are apart by 810 kilobases:

    >>> exons = BLAST.wellSeparatedExons(exons) # Keep exons separated by > 810KB
    Identifying exons separated by 810000 bases ...
    There are 575 exons

7. Finally we can use a function to save the obtained exons while making sure they are in frame. We need to use as additional arguments the genome file and output filename:

    >>> BLAST.storeExonsInFrame(exons, "pulled_seqs.fasta", "Bombyx_exons.fasta") 
    Storing exons ...
    A total of 575 exons are kept
    These exons have been stored in the file: Bombyx_exons.fasta


----------------------------
Validation of exon structure
----------------------------

We have now 575 single copy exons extracted from the *Bombyx mori* genome. Let's find
out whether these exons are conserved in other Arthropoda species.

For example we can compare these 575 exons with the genome of the monarch butterfly
*Danaus plexippus*.

^^^^^^^^^^^^^^^^^^
*Danaus plexippus*
^^^^^^^^^^^^^^^^^^

1. Download the version two of the monarch butterfly genome from here: http://danaus.genomeprojectsolutions-databases.com/Genome_seq_stats.html
2. Extract the genome as FASTA file using ``gunzip``:

   * ``gunzip Dp_genome_v2.fasta.gz``

3. Do a blastn of our Long Exons against the *Danaus* genome:

    >>> BLAST.blastn("Bombyx_exons.fasta", "Dp_genome_v2.fasta");
    ...
    BLASTn finished!
    The BLAST results were written in to the file Bombyx_exons_blastn_out.csv
    
4. We need to parse the output blast table and extract the exons from *Danaus* that are longer than 300bp and are homologous to the exons of *Bombyx mori*.

    >>> BLAST.blastParser("Bombyx_exons_blastn_out.csv", "Dp_genome_v2.fasta", "Danaus_exons.fasta", sp_name="Danaus")
    Reading files ...
    Parsing BLAST table ...
    A total of 158 sequences passed the thresholds.
    They have been stored in the file: Danaus_exons.fasta

The parameter ``sp_name`` is important as it will be used as part of the exons IDs.
 

^^^^^^^^^^^^^^^^^^^^^^
*Heliconius melpomene*
^^^^^^^^^^^^^^^^^^^^^^

1. We can continue finding homologous exons in other related butterflies. For example *Heliconius melpomene*.
2. Download the genome from here: http://metazoa.ensembl.org/Heliconius_melpomene/Info/Index
3. Extract the genome as FASTA file:

    * ``gunzip Heliconius_melpomene.Hmel1.17.dna_rm.toplevel.fa.gz``
    * ``mv Heliconius_melpomene.Hmel1.17.dna_rm.toplevel.fa Heliconius_genome.fa``

4. BLASTn the *Bombyx mori* exons against the *Heliconius* genome:

    >>> BLAST.blastn("Bombyx_exons.fasta", "Heliconius_genome.fa");
    ...
    BLASTn finished!
    The BLAST results were written in to the file  Bombyx_exons_blastn_out.csv

5. Parse the blast table, extract the exon sequences and save them to a file:

    >>> BLAST.blastParser("Bombyx_exons_blastn_out.csv", "Heliconius_genome.fa", "Heliconius_exons.fasta", sp_name="Heliconius")
    Reading files ...
    Parsing BLAST table ...
    A total of 145 sequences passed the thresholds.
    They have been stored in the file: Heliconius_exons.fasta
    
	    
^^^^^^^^^^^^^^^
*Manduca sexta*
^^^^^^^^^^^^^^^

1. Repeating the procedure for the *tobacco hornworm*.
2. Download the genome from ftp://ftp.bioinformatics.ksu.edu/pub/Manduca/
3. We downloaded the file ``Msex05162011.genome.fa``.
4. Blasted the *Bombyx mori* exons against the *Manduca* genome:

    >>> BLAST.blastn("Bombyx_exons.fasta", "Msex05162011.genome.fa")
    ...
    BLASTn finished!
    The BLAST results were written in to the file  Bombyx_exons_blastn_out.csv
5. Parsing the output blast table:

    >>> BLAST.blastParser("Bombyx_exons_blastn_out.csv", "Msex05162011.genome.fa", "Manduca_exons.fasta", sp_name="Manduca")
    Reading files ...
    Parsing BLAST table ...
    A total of 219 sequences passed the thresholds.
    They have been stored in the file: Manduca_exons.fasta
    

-----------
Small break
-----------

A **quick summary** of the work so far:

#. We obtained a list of orthologous and single copy genes by parsing the dataset for Arthopod genes from OrthoDB_.
#. From those genes, we took the exon sequences for *Bombyx mori* from its Coding DNA Sequences (CDS) from silkdb.org_.
#. We want to be sure that there are no introns inside our candidate exons. So we blasted the CDS sequences against the *Bombyx mori* genome.
#. We filtered those exons that were longer than 300 bp, were separated by 810 kilobases and got them inframe.
#. We did massive blasting of these selection of exons against genomes of related species: *Danaus plexippus*, *Heliconius melpomene* and *Manduca sexta*.
#. We got one FASTA file with the homologous regions for each species genome.
#. Now, we will proceed to align all those homologous exons in order to design primers.
#. Thus, we will be albe to sequence these exons accross a wide range of species in the order Lepidoptera.

--------------
Exon Alignment
--------------

We will use our module ``MUSCLE`` to do the alignment. We need to use as input a python list of the filenames that contain the exons of each species. 
All aligned sequences will be written into a folder called ``alignments`` as FASTA files (one file per exon).

.. warning:: 
  In the list of files, we will put **FIRST** the file for *Bombyx mori*, so that it will be used as "master" file. This is because the script will look for sequences in other files that appear in the file for *Bombyx mori*.

Example:

    >>> from pyphylogenomics import MUSCLE
    >>> files = ['Bombyx_exons.fasta', 'Danaus_exons.fasta','Heliconius_exons.fasta','Manduca_exons.fasta']
    >>> MUSCLE.batchAlignment(files)
    ...
    Pooling gene BGIBMGA000851:1-597
    Pooling gene BGIBMGA010204:1-516
    132 alignments have been saved in the folder "alignments"

.. warning::
    It is always recommended to check by eye every alignment that has been generated by any software. Once you are sure that the alignment is correct, we can continue with the analysis.


-------------
Primer design
-------------
Now that we have our exons/genes from several species (*Bombyx*, *Manduca*, *Danaus* and *Heliconius* in this example), we can design primers in order to sequence these genes across a wide range of butterflies and/or moths.

Since we have 132 candidate genes to design primers for, we can automate the primer design using a nice tool available in **PyPhyloGenomics**.

The function ``designPrimers`` will send an alignment to primers4clades_ along (with some parameters) and do a request for primer design. This function will return the degenerate primers as estimated by primers4clades_.

It is recommended that you enter your email as one of the parameters so that primers4clades_ can send you an email with very detailed results for you to insect. Just in case you don't provide your email, the very detailed results will be saved in the same folder of your alignments as HTML files.

.. warning:: Please keep in mind that the ``designPrimers`` function will return very little data, i.e. only the forward and reverse primers for an alignment. But it might be necessary that you inspect the detailed information saved as HTML files or the emails sent to you by primers4clades_. You will receive one email for each alignment.

Before primer design it could be useful to trim off the ends of the sequences
so that all sequences will have the same length:

    >>> from pyphylogenomics import MUSCLE
    >>> MUSCLE.bluntSplicer("alignments/") # folder_path containing the FASTA file alignments

This will produce FASTA files ending in **"_bluntlySpliced.fasta"**. You may want to remove the old 
unspliced FASTA files before doing primer design.

Automated primer design via primers4clades_:

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
    ...
    Done.
    All primers have been saved in the file "primers.fasta"
    
All primers will be saved to a file (``primers.fasta``). However, it is recommended that you study the very detailed
output saved into your ``alignments`` folder as HTML files so that you can decide to use these primers or not.

.. _primers4clades: http://floresta.eead.csic.es/primers4clades/#0
