#PyPhyloGenomics
A package to work on Phylogenomics.

[[In development.]]

###Developers
* Carlos Peña (email: carlos.pena@utu.fi)
* Victor Solis
* Pavel Matos
* Chris Wheat

###Install dependencies:

####Parallel Python (pp):
* If you are using Ubuntu Linux or related:

	sudo apt-get install python-pp

* Otherwise, install the distribution of `pp` included in this package:

	unzip pp-1.6.4.zip
	cd pp-1.6.4
	python setup.py install


###Installing PyPhyloGenomics
Use `setup.py`:

    python setup.py build  
    python setup.py install

###Reading the documentation
After installling:

    cd doc  
    make html

Then open the file `_build/html/index.html` in your web-browser.


