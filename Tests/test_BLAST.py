import unittest
import os
from pyphylogenomics import BLAST
from pyphylogenomics import OrthoDB


class BLASTTest(unittest.TestCase):

    def setUp(self):
        self.genes = OrthoDB.single_copy_genes("OrthoDB/OrthoDB6_Arthropoda_tabtext.csv", \
                    "Bombyx mori")

    def test_get_cds(self):
        BLAST.get_cds(self.genes, "BLAST/silkcds.fa")
        f = open("pulled_seqs.fasta", "r")
        result = len(f.read())
        f.close()
        """Extracting genes and saving them as fasta file"""
        self.assertEqual(result, 83216)
        os.remove("pulled_seqs.fasta")





if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
