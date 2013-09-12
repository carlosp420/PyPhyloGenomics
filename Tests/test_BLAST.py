import unittest
import os
from pyphylogenomics import BLAST
from pyphylogenomics import OrthoDB


class BLASTTest(unittest.TestCase):

    def setUp(self):
        self.genes = OrthoDB.single_copy_genes("OrthoDB/OrthoDB6_Arthropoda_tabtext.csv", \
                    "Bombyx mori")
        self.genome = "BLAST/silkcds.fa"

    def test_get_cds(self):
        BLAST.get_cds(self.genes, "BLAST/silkcds.fa")
        f = open("pulled_seqs.fasta", "r")
        result = len(f.read())
        f.close()
        """Extracting genes and saving them as fasta file"""
        self.assertEqual(result, 83216)
        os.remove("pulled_seqs.fasta")

    def test_makeblastdb_true(self):
        mask = True
        BLAST.makeblastdb(self.genome, mask)
        for name in os.listdir("BLAST/"):
            if name[:10] == "silkcds.fa" and len(name) > 10:
                os.remove("BLAST/" + name);
                result = "true"
        self.assertEqual(result, "true");

    def test_makeblastdb_false(self):
        mask = False
        BLAST.makeblastdb(self.genome, mask)
        for name in os.listdir("BLAST/"):
            if name[:10] == "silkcds.fa" and len(name) > 10:
                os.remove("BLAST/" + name);
                result = "true"
        self.assertEqual(result, "true");

    def test_blastn(self)
        BLAST.blastn("BLAST/query.fas", "BLAST/silkcds.fa");
        




if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
