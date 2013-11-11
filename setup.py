#from distutils.core import setup
from distribute_setup import use_setuptools
use_setuptools()
from setuptools import setup, find_packages


setup(
		name='PyPhyloGenomics', 
		version='0.3.10',
		packages = ['pyphylogenomics'],
		author="Carlos Pena, Victor Solis, Pavel Matos, Chris Wheat",
		author_email="mycalesis@gmail.com",
		maintainer="Carlos Pena, Victor Solis, Pavel Matos, Chris Wheat",
		maintainer_email="mycalesis@gmail.com",
		contact="Carlos Pena",
		contact_email="mycalesis@gmail.com",
		license="GPL v3",
		description="Tools to work in phylogenomics, from NSG group http://nymphalidae.utu.fi",
		long_description=open('README.txt').read(),
		keywords="DNA, genomics, genomes, phylogenetics, genes",
		url="http://carlosp420.github.com/PyPhyloGenomics/",
		download_url = "http://carlosp420.github.com/PyPhyloGenomics/",
		classifiers=[
					"Programming Language :: Python",
					("Topic :: Scientific/Engineering :: Bio-Informatics"),
					("Intended Audience :: Science/Research"),
					("License :: OSI Approved :: GNU General Public License v2 (GPLv2)"),
					("Operating System :: OS Independent"),
					("Environment :: Console"),
					],

		install_requires=['biopython','pp','requests','beautifulsoup4']
		);
