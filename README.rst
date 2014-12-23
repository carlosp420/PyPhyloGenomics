|Build Status| |Coveralls|

PyPhyloGenomics A package to work on Phylogenomics
==================================================

Developers
----------

-  Carlos Peña (email: mycalesis@gmail.com)
-  Victor Solis
-  Pavel Matos
-  Chris Wheat

Installing PyPhyloGenomics
--------------------------

PyPhyloGenomics has been developed in Python v2.7. The installer of
PyPhyloGenomics will try to download and install all its dependencies as
well.
Download the latest version from here:
https://github.com/carlosp420/PyPhyloGenomics/releases

To install PyPhyloGenomics use ``setup.py``:

::

    python setup.py build  
    python setup.py install

If it fails you can install the dependencies manually:

Install dependencies:
---------------------

requests:
~~~~~~~~~

The package ``requests`` from
`here <http://docs.python-requests.org/en/latest/user/install/>`__. Or
try:

::

    sudo apt-get install python-requests

BioPython:
~~~~~~~~~~

Download and install from `here <http://biopython.org/wiki/Download>`__.
Or:

::

    sudo apt-get install python-biopython

BeatutifulSoup
~~~~~~~~~~~~~~

Download and install from
`here <http://www.crummy.com/software/BeautifulSoup/>`__. Or:

::

    sudo apt-get install python-bs4

MUSCLE
~~~~~~

It is necessary that you install MUSCLE so that PyPhyloGenomics can use
it to align sequences. Download and install from
`here <http://www.drive5.com/muscle/downloads.htm>`__.

If you are using Windows you can download the executable file
**muscle3.8.31\_i86win32.exe** and save it in your Python folder (C:27)
as **muscle.exe**.

BLAST
~~~~~

Download and install the BLAST+ executables from the `NCBI
website <http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download>`__.
Or try:

::

    sudo apt-get install ncbi-blast+

fastx-toolkit
~~~~~~~~~~~~~

Download and install from
`here <http://hannonlab.cshl.edu/fastx_toolkit/>`__. Or:

::

    sudo apt-get install fastx-toolkit

Reading PyPhyloGenomics' documentation:
---------------------------------------

Read the online documentation here:
http://pyphylogenomics.readthedocs.org/

Or, after installling do the following:

::

    cd doc  
    make html

Then open the file ``_build/html/index.html`` in your web-browser.

Reproduce our analysis:
-----------------------

You can reproduce all the steps detailed in our [Getting started with
PyPhylogenomics] guide. Just use the command line in the same folder as
the ``Makefile`` and type make ???

.. |Build Status| image:: https://travis-ci.org/carlosp420/PyPhyloGenomics.png?branch=master
   :target: https://travis-ci.org/carlosp420/PyPhyloGenomics

.. |Coveralls| image:: https://coveralls.io/repos/carlosp420/PyPhyloGenomics/badge.png?branch=master
  :target: https://coveralls.io/r/carlosp420/PyPhyloGenomics?branch=master
