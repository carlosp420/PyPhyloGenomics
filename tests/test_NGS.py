import unittest
import os
import os.path
from pyphylogenomics import NGS
import subprocess
import glob
from Bio import SeqIO
import shutil


class NGSTest(unittest.TestCase):

    def setUp(self):
        self.cwd = os.path.dirname(__file__)

    def test_prepare_data(self):
        ionfile = os.path.join(self.cwd, "NGS", "ion_file.fastq")
        NGS.prepare_data(ionfile, 8)

        expected_file1 = os.path.join("data", "modified", "wrk_ionfile.fasta")
        expected_file2 = os.path.join("data", "modified", "wrk_ionfile.fastq")
        self.assertTrue(os.path.isfile(expected_file1))
        self.assertTrue(os.path.isfile(expected_file2))
        os.remove(expected_file1)
        os.remove(expected_file2)

    def test_filter_reads(self):
        folder = os.path.join(self.cwd, "NGS")

        for i in glob.glob(os.path.join("NGS", "gene*")):
            os.remove(i)

        ion_chunk = os.path.join(self.cwd, "NGS", "reaa.fastq")
        blast_chunk = os.path.join(self.cwd, "NGS", "reaa.csv")
        NGS.filter_reads(ion_chunk, blast_chunk, folder)

        cmd = "cat " + os.path.join(self.cwd, "NGS", "gene*")
        cmd += " | grep -c '^@[0-9]\+'"
        p = subprocess.check_output(cmd, shell=True)

        for i in glob.glob(os.path.join(self.cwd, "NGS", "gene*")):
            os.remove(i)
        self.assertEqual(int(p.strip()), 23)

    def test_parse_blast_results(self):
        # It should work using fasta files
        blast_table = os.path.join(self.cwd, "NGS", "blast_table.csv")
        ion_file = os.path.join(self.cwd, "NGS", "ion_file.fastq")

        NGS.parse_blast_results(blast_table, ion_file)
        result = glob.glob("output/gene*")
        self.assertEqual(len(result), 21)
        shutil.rmtree("output")

    def test_split_ionfile_by_results(self):
        ion_file = os.path.join(self.cwd, "NGS/ion_file.fastq")
        blast_chunk = os.path.join(self.cwd, "NGS/_reaa.csv")

        NGS.split_ionfile_by_results(ion_file, blast_chunk)
        cmd = "grep -c '^@' " + os.path.join(self.cwd, "NGS", "_reaa.fastq")
        p = subprocess.check_output(cmd, shell=True)
        self.assertEqual(int(p.strip()), 1001)
        os.remove(os.path.join(self.cwd, "NGS", "_reaa.fastq"))

    def test_prune(self):
        folder = os.path.join(self.cwd, "NGS")

        blast_data = []
        f = open(os.path.join(folder, "blast_data.csv"), "r")
        tmp = f.readlines()
        f.close()
        for i in tmp:
            blast_data.append(i.strip())

        seq_record = SeqIO.parse(os.path.join(folder, "seq_record.fastq"), "fastq")

        ion_id = "3856"
        min_aln_length = "40"

        result = NGS.prune(folder, blast_data, seq_record, ion_id,
                           min_aln_length)
        # It should drop on seq_record from the blast_data
        self.assertEqual(len(result), 998)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
