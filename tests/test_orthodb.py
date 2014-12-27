import os
import unittest

from pyphylogenomics import OrthoDB


class OrthoDBTest(unittest.TestCase):
    def setUp(self):
        self.cwd = os.path.dirname(__file__)
        self.in_file = os.path.join(self.cwd, "OrthoDB", "OrthoDB6_Arthropoda_tabtext.csv")

    def test_single_copy_genes(self):
        """We are using only part of the original OrthoDB6 file"""
        result = OrthoDB(self.in_file, 'Bombyx mori')
        genes = result.single_copy_genes
        self.assertEqual(len(genes), 325)

    def test_copies_per_gene(self):
        result = OrthoDB.copies_per_gene(self.in_file)
        number_copies = result['Bombyx mori', 'BGIBMGA000894', 'EOG600001']
        self.assertEqual(2, number_copies)


    def test_single_copy_in_species(self):
        result = OrthoDB.single_copy_in_species(
            os.path.join(self.cwd, "OrthoDB", "OrthoDB6_Arthropoda_tabtext.csv"), "Znev_18776")
        self.assertEqual(result[0], "Zootermopsis nevadensis")


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
