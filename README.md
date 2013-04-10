#PyPhyloGenomics
A package to work on Phylogenomics.

[[In development.]]

###Developers
* Carlos Peña (email: carlos.pena@utu.fi)
* Victor Solis
* Pavel Matos
* Chris Wheat


###Installing PyPhyloGenomics
The installer of PyPhyloGenomics will try to download and install all its dependencies as well. 
To install PyPhyloGenomics use `setup.py`:

    python setup.py build  
    python setup.py install

If it fails you can install the dependencies manually:

###Install dependencies:

####Parallel Python (pp):
If you are using Ubuntu Linux or related:

    sudo apt-get install python-pp

Otherwise, [download](http://www.parallelpython.com/content/view/15/30/) the source code and install `pp`:

    unzip pp-1.6.4.zip
    cd pp-1.6.4
    python setup.py install

####BioPython:
Download and install from [here](http://biopython.org/wiki/Download).

###Reading PyPhyloGenomics' documentation:
After installling:

    cd doc  
    make html

Then open the file `_build/html/index.html` in your web-browser.


