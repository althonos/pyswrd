import unittest

from .. import HeuristicFilter, Sequences, search
from . import data


class TestSearch(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        cls.o74807 = Sequences()
        for record in data.load_records("O74807.fasta"):
            cls.o74807.append(record.seq)
        cls.sprot15 = Sequences()
        cls.names = []
        for record in data.load_records("uniprot_sprot15.fasta"):
            cls.sprot15.append(record.seq)
            cls.names.append(record.id)
        cls.sprot196 = Sequences()
        for record in data.load_records("uniprot_sprot196.fasta"):
            cls.sprot196.append(record.seq)
        cls.sprot12071 = Sequences()
        for record in data.load_records("uniprot_sprot12071.fasta"):
            cls.sprot12071.append(record.seq)

    def test_search(self):
        # Results from SWORD:
        # Query id,Subject id,% identity,alignment length,mismatches,gap openings,q. start,q. end,s. start,s. end,e-value,score
        # sp|O74807|YGNG_SCHPO    sp|Q6GZW6|009L_FRG3G    23      30      23      0       67      96      415     444     1.32    35
        # sp|O74807|YGNG_SCHPO    sp|Q91G88|006L_IIV6     29      41      23      2       2       39      131     168     1.64    34
        # sp|O74807|YGNG_SCHPO    sp|Q6GZW8|007R_FRG3G    40      15      9       0       7       21      29      43      2.42    32
        # sp|O74807|YGNG_SCHPO    sp|Q6GZX4|001R_FRG3G    31      36      24      1       59      94      133     167     2.76    32
        # sp|O74807|YGNG_SCHPO    sp|Q6GZX3|002L_FRG3G    24      82      43      5       15      85      103     176     4.92    30
        # sp|O74807|YGNG_SCHPO    sp|Q197F5|005L_IIV3     40      15      9       0       6       20      114     128     6.23    29
        # sp|O74807|YGNG_SCHPO    sp|Q197F2|008L_IIV3     29      41      27      1       4       42      261     301     6.54    29

        hits = list(search(self.o74807, self.sprot15))
        self.assertEqual(len(hits), 7)
        hits.sort(key=lambda h: h.evalue)

        self.assertAlmostEqual(hits[0].evalue, 1.32, places=2)
        self.assertAlmostEqual(hits[1].evalue, 1.64, places=2)
        self.assertAlmostEqual(hits[2].evalue, 2.42, places=2)
        self.assertAlmostEqual(hits[3].evalue, 2.76, places=2)
        self.assertAlmostEqual(hits[4].evalue, 4.92, places=2)
        self.assertAlmostEqual(hits[5].evalue, 6.23, places=2)
        self.assertAlmostEqual(hits[6].evalue, 6.54, places=2)

        self.assertEqual(hits[0].score, 35)
        self.assertEqual(hits[1].score, 34)
        self.assertEqual(hits[2].score, 32)
        self.assertEqual(hits[3].score, 32)
        self.assertEqual(hits[4].score, 30)
        self.assertEqual(hits[5].score, 29)
        self.assertEqual(hits[6].score, 29)

        self.assertEqual(hits[0].target_index, 13)
        self.assertEqual(hits[1].target_index, 8)
        self.assertEqual(hits[2].target_index, 10)
        self.assertEqual(hits[3].target_index, 0)
        self.assertEqual(hits[4].target_index, 1)
        self.assertEqual(hits[5].target_index, 6)
        self.assertEqual(hits[6].target_index, 12)
        
        # import pyopal
        # hits = pyopal.align( self.o74807[0], list(self.sprot15), score_matrix=pyopal.ScoreMatrix.aa("BLOSUM62"), gap_extend=1, gap_open=10 )

        # for hit in hits:

        #     query = self.o74807[0]
        #     target = self.sprot15[hit.target_index]
        #     # print(self.names[hit.target_index])

        #     # print(query)
        #     # print(target)

        #     # print(hit.target_index, hit.evalue, hit.score, sep="\t")
        #     print(hit.target_index, hit.score, sep="\t")



    # def test_search(self):
    #     hits = list(search(self.o74807, self.sprot15, gap_open=10, gap_extend=1, scorer_name="BLOSUM50", algorithm="sw"))
    #     self.assertEqual(len(hits), 10)
        
    #     hits.sort(key=lambda h: h.evalue)

    #     self.assertAlmostEqual(hits[0].evalue, 1.70e-5)
    #     self.assertAlmostEqual(hits[1].evalue, 1.87e-5)
    #     self.assertAlmostEqual(hits[2].evalue, 2.34e-5)

    #     for hit in hits:

    #         # query = self.o74807[0]
    #         # target = self.sprot15[hit.target_index]
    #         # print(self.names[hit.target_index])

    #         # print(query)
    #         # print(target)

    #         # print(hit.result.query_start, hit.result.query_end, hit.result.target_start, hit.result.target_end, hit.evalue, hit.score, sep="\t")


