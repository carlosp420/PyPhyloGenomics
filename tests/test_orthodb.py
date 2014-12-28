import os
import unittest

from pyphylogenomics import OrthoDB


class OrthoDBTest(unittest.TestCase):
    def setUp(self):
        self.cwd = os.path.dirname(__file__)
        self.in_file = os.path.join(self.cwd, "OrthoDB", "OrthoDB6_Arthropoda_tabtext.csv")
        self.my_orthodb = OrthoDB(self.in_file, 'Bombyx mori')
        self.outfile = os.path.join(self.cwd, "OrthoDB", "output.txt")

    def test_single_copy_genes(self):
        """We are using only part of the original OrthoDB6 file"""
        genes = self.my_orthodb.single_copy_genes
        self.assertEqual(len(genes), 325)

    def test_copies_per_gene(self):
        result = self.my_orthodb._copies_per_gene()
        number_copies = result['Bombyx mori', 'BGIBMGA000894', 'EOG600001']
        self.assertEqual(2, number_copies)

    def test_single_copy_in_species(self):
        result = self.my_orthodb._single_copy_in_species("Znev_18776")
        self.assertEqual(result[0], "Zootermopsis nevadensis")

    def test_copies_per_gene_table(self):
        self.my_orthodb._single_copy_in_species("Znev_18776")
        self.my_orthodb._copies_per_gene_table(self.outfile)
        with open(self.outfile, "r") as handle:
            self.assertTrue(handle.read().startswith('Genes\tZootermopsis nevadensis'))

    def tearDown(self):
        if os.path.isfile(self.outfile):
            os.remove(self.outfile)


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
