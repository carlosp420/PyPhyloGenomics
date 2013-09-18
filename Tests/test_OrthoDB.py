import unittest
from pyphylogenomics import OrthoDB


class OrthoDBTest(unittest.TestCase):

    def test_single_copy_genes(self):
        """We are using only part of the original OrthoDB6 file"""
        result = OrthoDB.single_copy_genes("OrthoDB/OrthoDB6_Arthropoda_tabtext.csv", \
                    "Bombyx mori")
        self.assertEqual(len(result), 446)

    def test_single_copy_in_species(self):
        result = \
        OrthoDB.single_copy_in_species("OrthoDB/OrthoDB6_Arthropoda_tabtext.csv", \
                    "Znev_18776")
        self.assertEqual(result[0], "Zootermopsis nevadensis")


if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
