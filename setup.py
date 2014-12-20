#!/usr/bin/env python
# -*- coding: utf-8 -*-


try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


readme = open('README.rst').read()
history = open('HISTORY.rst').read().replace('.. :changelog:', '')

requirements = [
    # TODO: put package requirements here
]

test_requirements = [
    # TODO: put package test requirements here
]

setup(
    name='pyphylogenomics',
    version='0.3.12',
    description='Tools to work in phylogenomics, from NSG group http://nymphalidae.utu.fi',
    long_description=readme + '\n\n' + history,
    author='Carlos Pena, Victor Solis, Pavel Matos, Chris Wheat',
    author_email='mycalesis@gmail.com',
		maintainer="Carlos Pena, Victor Solis, Pavel Matos, Chris Wheat",
		maintainer_email="mycalesis@gmail.com",
		contact="Carlos Pena",
		contact_email="mycalesis@gmail.com",
    url='https://github.com/carlosp420/PyPhyloGenomics',
    packages=[
        'pyphylogenomics',
    ],
    package_dir={'pyphylogenomics':
                 'pyphylogenomics'},
    include_package_data=True,
    install_requires=requirements,
    license="GPL v3",
    zip_safe=False,
    keywords='DNA, genomics, genomes, phylogenetics, genes',
    classifiers=[
					"Programming Language :: Python",
					("Topic :: Scientific/Engineering :: Bio-Informatics"),
					("Intended Audience :: Science/Research"),
					("License :: OSI Approved :: GNU General Public License v2 (GPLv2)"),
					("Operating System :: OS Independent"),
					("Environment :: Console"),
					],
    test_suite='tests',
    tests_require=test_requirements
)
