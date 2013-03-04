from distutils.core import setup

setup(
		name='PyPhyloGenomics', 
		version='0.2.0',
		packages = ['pyphylogenomics'],
		author="Carlos Pena, Victor Solis, Pavel Matos, Chris Wheat",
		author_email="mycalesis@gmail.com",
		maintainer="Carlos Pena, Victor Solis, Pavel Matos, Chris Wheat",
		maintainer_email="mycalesis@gmail.com",
		contact="Carlos Pena",
		contact_email="mycalesis@gmail.com",
		license="GPL v2",
		description="Tools to work in phylogenomics, from NSG group http://nymphalidae.utu.fi",
		long_description=open('README.txt').read(),
		keywords="DNA, genomics, genomes, phylogenetics, genes",
		url="http://nymphalidae.utu.fi",
		download_url = "https://github.com/carlosp420/PyPhyloGenomics/",
		classifiers=[
					"Programming Language :: Python",
					("Topic :: Scientific/Engineering :: Bio-Informatics"),
					("Intended Audience :: Science/Research"),
					("License :: OSI Approved :: GNU General Public License v2 (GPLv2)"),
					("Operating System :: OS Independent"),
					("Environment :: Console"),
					],

		requires=['biopython','MySQLdb']
		);
