from distutils.core import setup

setup(name='PyPhyloGenomics', version='0.1.0',
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
		classifiers=[
					"Programming Language :: Python",
					("Topic :: Scientific/Engineering :: Bio-Informatics"),
					("Intended Audience :: Science/Research"),
					("Environment :: Console"),
					],

		install_requires=['biopython','MySQLdb']);
