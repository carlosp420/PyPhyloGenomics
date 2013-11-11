import unittest
import os
import os.path
from pyphylogenomics import NGS
import subprocess
import glob


class NGSTest(unittest.TestCase):

    def setUp(self):
        pass

    def test_filter_reads(self):
        folder = "NGS"

        for i in glob.glob(os.path.join("NGS", "gene*")):
            os.remove(i)

        ion_chunk = os.path.join("NGS", "_reaa.fastq")
        blast_chunk = os.path.join("NGS", "_reaa.csv")
        NGS.filter_reads(ion_chunk, blast_chunk, folder)

        cmd = "cat " + os.path.join("NGS", "gene*")
        cmd += " | grep -c '^@'"
        p = subprocess.check_output(cmd, shell=True)

        for i in glob.glob(os.path.join("NGS", "gene*")):
            os.remove(i)
        self.assertEqual(int(p.strip()), 23)

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)

