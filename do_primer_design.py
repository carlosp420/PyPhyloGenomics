#!/usr/bin/env python
from pyphylogenomics import MUSCLE

"""
Automated primer design via primers4clades:
"""
folder = "alignments"   # folder containing the FASTA file alignments
tm = "55"               # annealing temperature
min_amplength = "250"   # minimium amplicon length
max_amplength = "500"   # maximum amplicon length
gencode = "universal"   # see below for all available genetic codes
mode  = "primers"
clustype = "dna"
amptype = "dna_GTRG"    # substitution model used to estimate phylogenetic information
email = "mycalesis@gmail.com"   # primer4clades will send you an email with very detailed results

MUSCLE.designPrimers(folder, tm, min_amplength, max_amplength, gencode, mode, clustype, amptype, email)

