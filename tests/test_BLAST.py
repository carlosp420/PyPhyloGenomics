import subprocess
import unittest
import os

from pyphylogenomics import BLAST
from pyphylogenomics import OrthoDB


class BLASTTest(unittest.TestCase):

    def setUp(self):
        self.cwd = os.path.dirname(__file__)
        self.genes = OrthoDB.single_copy_genes(
            os.path.join(self.cwd, "OrthoDB", "OrthoDB6_Arthropoda_tabtext.csv"),
            "Bombyx mori"
        )
        self.genome = os.path.join(self.cwd, "BLAST", "silkcds.fa")

    def test_get_cds(self):
        """Extracting genes and saving them as fasta file"""
        BLAST.get_cds(self.genes, self.genome)
        f = open("pulled_seqs.fasta", "r")
        result = len(f.read())
        f.close()
        self.assertEqual(result, 63178)
        os.remove("pulled_seqs.fasta")

    def test_makeblastdb_true(self):
        mask = True
        BLAST.makeblastdb(self.genome, mask)
        for name in os.listdir(os.path.join(self.cwd, "BLAST") + "/"):
            if name[:10] == "silkcds.fa" and len(name) > 10:
                os.remove(self.cwd + "/BLAST/" + name)
                result = "true"
        self.assertEqual(result, "true")

    def test_makeblastdb_false(self):
        mask = False
        BLAST.makeblastdb(self.genome, mask)
        for name in os.listdir(self.cwd + "/BLAST/"):
            if name[:10] == "silkcds.fa" and len(name) > 10:
                os.remove(self.cwd + "/BLAST/" + name)
                result = "true"
        self.assertEqual(result, "true")

    def test_do_blast(self):
        f = os.path.join(self.cwd, 'BLAST', 'query.fas')
        genome = os.path.join(self.cwd, 'BLAST', 'silkcds.fa')
        BLAST.makeblastdb(genome)

        command = 'blastn -query ' + f + ' -db ' + genome + ' -task blastn ' \
                  + '-evalue 0.0001 ' + ' -out ' + f + '_out.csv' + ' -num_threads 1 ' \
                  + ' -outfmt 10'
        BLAST.do_blast(command)
        output = open(os.path.join(self.cwd, 'BLAST', 'query.fas_out.csv'), "r")
        for line in output:
            self.assertTrue('BGIBMGA000001' in line)
            break

        for name in os.listdir(self.cwd + "/BLAST/"):
            if name[:10] == "silkcds.fa" and len(name) > 10:
                os.remove(self.cwd + "/BLAST/" + name)
        os.remove(os.path.join(self.cwd, 'BLAST', 'query.fas_out.csv'))

    def test_blastn(self):
        BLAST.blastn(self.cwd + "/BLAST/query.fas", self.cwd + "/BLAST/silkcds.fa")
        file = open(self.cwd + "/BLAST/query_blastn_out.csv", "r")
        for line in file:
            result = line.split(",")[1]
            break
        for name in os.listdir(self.cwd + "/BLAST/"):
            if name[:10] == "silkcds.fa" and len(name) > 10:
                os.remove(self.cwd + "/BLAST/" + name)
        self.assertEqual(result, "BGIBMGA000001-TA")

    def test_blastn_big_query_file(self):
        query_file = os.path.join(self.cwd, "BLAST", "query_big.fas.gz")

        cmd = "gunzip " + query_file
        p = subprocess.check_call(cmd, shell=True)

        if p == 0:
            gunzipped_query_file = os.path.join(self.cwd, "BLAST", "query_big.fas")
            BLAST.blastn(gunzipped_query_file, self.cwd + "/BLAST/silkcds.fa")
            file = open(self.cwd + "/BLAST/query_big_blastn_out.csv", "r")
            for line in file:
                result = line.split(",")[1]
                break
            for name in os.listdir(self.cwd + "/BLAST/"):
                if name[:10] == "silkcds.fa" and len(name) > 10:
                    os.remove(self.cwd + "/BLAST/" + name)
            self.assertEqual(result, "BGIBMGA000001-TA")

            os.remove(os.path.join(self.cwd, "BLAST", "query_big_blastn_out.csv"))
            cmd = "gzip " + gunzipped_query_file
            p = subprocess.check_call(cmd, shell=True)
        else:
            raise Exception('test failed.')

    def test_getLargestExon(self):
        exons = BLAST.getLargestExon(
            self.cwd + "/BLAST/query_blastn_output1.csv",
            E_value=0.001, ident=98, exon_len=300,
        )
        result = len(exons)
        self.assertEqual(result, 36)

    def test_getLargestExon_output_has_headers(self):
        exons = BLAST.getLargestExon(
            self.cwd + "/BLAST/query_blastn_output2.csv",
            E_value=0.001, ident=98, exon_len=300,
        )
        result = len(exons)
        self.assertEqual(result, 38)

    def test_eraseFalsePosi(self):
        exons = BLAST.getLargestExon(
            self.cwd + "/BLAST/query_blastn_out.csv",
            E_value=0.001, ident=98, exon_len=300,
        )
        exons = BLAST.eraseFalsePosi(exons)
        result = len(exons)
        self.assertEqual(result, 3)

    def test_wellSeparatedExons(self):
        exons = BLAST.getLargestExon(
            os.path.join(self.cwd, "BLAST", "query_blastn_output3.csv"),
            E_value=0.001,
            ident=98,
            exon_len=300,
        )
        exons = BLAST.wellSeparatedExons(exons)
        result = len(exons)
        self.assertEqual(result, 3)

    def test_wellSeparatedExons_drop_one_exon(self):
        exons = BLAST.getLargestExon(
            os.path.join(self.cwd, "BLAST", "query_blastn_output4.csv"),
            E_value=0.001,
            ident=98,
            exon_len=300,
        )
        exons = BLAST.wellSeparatedExons(exons)
        result = len(exons)
        self.assertEqual(result, 2)

    def test_filterByMinDist(self):
        # The gene2 is too close to other genes
        genes_loci = [
            ('gene1', 1, 350),
            ('gene2', 360, 670),
            ('gene3', 821001, 821351),
        ]
        result = BLAST.filterByMinDist(genes_loci, 810000)
        self.assertEqual(['gene2'], result)

        # The function is not affected by order of genes
        genes_loci = [
            ('gene3', 821001, 821351),
            ('gene1', 1, 350),
            ('gene2', 360, 670),
        ]
        result = BLAST.filterByMinDist(genes_loci, 810000)
        self.assertEqual(['gene2'], result)

        # These genes are well separated already
        genes_loci = [
            ('gene1', 1, 350),
            ('gene3', 821001, 821351),
        ]
        result = BLAST.filterByMinDist(genes_loci, 810000)
        self.assertEqual([], result)
    # todo BLAST.storeExonsInFrame

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
