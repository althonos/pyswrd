import unittest

from .. import HeuristicFilter, Sequences
from . import data


class _TestFilter:
    @classmethod
    def setUpClass(cls):
        cls.o74807 = Sequences()
        if data.exists("O74807.fasta"):
            for record in data.load_records("O74807.fasta"):
                cls.o74807.append(record.seq)
        cls.sprot15 = Sequences()
        if data.exists("uniprot_sprot15.fasta"):
            for record in data.load_records("uniprot_sprot15.fasta"):
                cls.sprot15.append(record.seq)
        cls.sprot196 = Sequences()
        if data.exists("uniprot_sprot196.fasta"):
            for record in data.load_records("uniprot_sprot196.fasta"):
                cls.sprot196.append(record.seq)
        cls.sprot12071 = Sequences()
        if data.exists("uniprot_sprot12071.fasta"):
            for record in data.load_records("uniprot_sprot12071.fasta"):
                cls.sprot12071.append(record.seq)


class TestHeuristicFilter(_TestFilter, unittest.TestCase):

    maxDiff = None

    @unittest.skipUnless(data.exists("O74807.fasta"), "missing data file")
    @unittest.skipUnless(data.exists("uniprot_sprot15.fasta"), "missing data file")
    def test_chunks(self):
        f1 = HeuristicFilter(self.o74807, threads=1)
        r1 = f1.score(self.sprot15).finish()

        f2 = HeuristicFilter(self.o74807, threads=1)
        f2.score(self.sprot15.extract(range(10)))
        f2.score(self.sprot15.extract(range(10, len(self.sprot15))))
        r2 = f2.finish()

        self.assertEqual(len(r1.entries), len(r2.entries))
        for e1, e2 in zip(r1.entries, r2.entries):
            self.assertEqual(e1, e2)        

    @unittest.skipUnless(data.exists("O74807.fasta"), "missing data file")
    @unittest.skipUnless(data.exists("uniprot_sprot15.fasta"), "missing data file")
    def test_chunks_threads(self):
        f1 = HeuristicFilter(self.o74807, threads=4)
        r1 = f1.score(self.sprot15).finish()

        f2 = HeuristicFilter(self.o74807, threads=4)
        f2.score(self.sprot15.extract(range(10)))
        f2.score(self.sprot15.extract(range(10, len(self.sprot15))))
        r2 = f2.finish()

        self.assertEqual(len(r1.entries), len(r2.entries))
        for e1, e2 in zip(r1.entries, r2.entries):
            self.assertEqual(e1, e2)     




class TestK3(_TestFilter, unittest.TestCase):

    @unittest.skipUnless(data.exists("uniprot_sprot15.fasta"), "missing data file")
    @unittest.skipUnless(data.exists("uniprot_sprot196.fasta"), "missing data file")
    def test_sprot15_sprot196(self):
        f = HeuristicFilter(self.sprot15)
        f.score(self.sprot196)
        result = f.finish()
        self.assertEqual(len(result.indices[0]), 185)
        self.assertEqual(len(result.indices[1]), 193)
        self.assertEqual(len(result.indices[2]), 192)
        self.assertEqual(len(result.indices[3]), 184)
        self.assertEqual(len(result.indices[4]), 190)
        self.assertEqual(len(result.indices[5]), 76)
        self.assertEqual(len(result.indices[6]), 190)
        self.assertEqual(len(result.indices[7]), 177)
        self.assertEqual(len(result.indices[8]), 191)
        self.assertEqual(len(result.indices[9]), 132)
        self.assertEqual(len(result.indices[10]), 158)
        self.assertEqual(len(result.indices[11]), 192)
        self.assertEqual(len(result.indices[12]), 187)
        self.assertEqual(len(result.indices[13]), 195)
        self.assertEqual(len(result.indices[14]), 164)

    @unittest.skipUnless(data.exists("O74807.fasta"), "missing data file")
    @unittest.skipUnless(data.exists("uniprot_sprot15.fasta"), "missing data file")
    def test_o74807_sprot15(self):
        f = HeuristicFilter(self.o74807)
        f.score(self.sprot15)
        result = f.finish()
        self.assertEqual(
            result.indices[0], 
            [0, 1, 2, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14]
        )

    @unittest.skipUnless(data.exists("O74807.fasta"), "missing data file")
    @unittest.skipUnless(data.exists("uniprot_sprot196.fasta"), "missing data file")
    def test_o74807_sprot196(self):
        f = HeuristicFilter(self.o74807)
        f.score(self.sprot196)
        result = f.finish()
        self.assertEqual(len(result.indices[0]), 140)
        self.assertEqual(
            result.indices[0][:10], 
            [0, 1, 2, 4, 5, 6, 7, 8, 10, 11]
        )
        self.assertEqual(
            result.indices[0][-10:],
            [181, 187, 188, 189, 190, 191, 192, 193, 194, 195]
        )

    @unittest.skipUnless(data.exists("O74807.fasta"), "missing data file")
    @unittest.skipUnless(data.exists("uniprot_sprot12071.fasta"), "missing data file")
    def test_o74807_sprot12071(self):
        f = HeuristicFilter(self.o74807)
        f.score(self.sprot12071)
        result = f.finish()
        self.assertEqual(len(result.indices[0]), 10088)
        self.assertEqual(
            result.indices[0][:10], 
            [0, 1, 2, 4, 5, 6, 7, 8, 10, 11]
        )
        self.assertEqual(
            result.indices[0][-10:],
            [12061, 12062, 12063, 12064, 12065, 12066, 12067, 12068, 12069, 12070]
        )


class TestK4(_TestFilter, unittest.TestCase):

    @unittest.skipUnless(data.exists("uniprot_sprot15.fasta"), "missing data file")
    @unittest.skipUnless(data.exists("uniprot_sprot196.fasta"), "missing data file")
    def test_sprot15_sprot196(self):
        f = HeuristicFilter(self.sprot15, kmer_length=4)
        f.score(self.sprot196)
        result = f.finish()
        self.assertEqual(len(result.indices[0]), 196)
        self.assertEqual(len(result.indices[1]), 196)
        self.assertEqual(len(result.indices[2]), 196)
        self.assertEqual(len(result.indices[3]), 195)
        self.assertEqual(len(result.indices[4]), 196)
        self.assertEqual(len(result.indices[5]), 177)
        self.assertEqual(len(result.indices[6]), 195)
        self.assertEqual(len(result.indices[7]), 196)
        self.assertEqual(len(result.indices[8]), 196)
        self.assertEqual(len(result.indices[9]), 183)
        self.assertEqual(len(result.indices[10]), 191)
        self.assertEqual(len(result.indices[11]), 196)
        self.assertEqual(len(result.indices[12]), 196)
        self.assertEqual(len(result.indices[13]), 196)
        self.assertEqual(len(result.indices[14]), 176)

    @unittest.skipUnless(data.exists("O74807.fasta"), "missing data file")
    @unittest.skipUnless(data.exists("uniprot_sprot15.fasta"), "missing data file")
    def test_o74807_sprot15(self):
        f = HeuristicFilter(self.o74807, kmer_length=4)
        f.score(self.sprot15)
        result = f.finish()
        self.assertEqual(
            result.indices[0],
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14],
        )

    @unittest.skipUnless(data.exists("O74807.fasta"), "missing data file")
    @unittest.skipUnless(data.exists("uniprot_sprot196.fasta"), "missing data file")
    def test_o74807_sprot196(self):
        f = HeuristicFilter(self.o74807, kmer_length=4)
        f.score(self.sprot196)
        result = f.finish()
        self.assertEqual(set(range(196)) - set(result.indices[0]), {193, 110, 159})

    @unittest.skipUnless(data.exists("O74807.fasta"), "missing data file")
    @unittest.skipUnless(data.exists("uniprot_sprot12071.fasta"), "missing data file")
    def test_o74807_sprot12071(self):
        f = HeuristicFilter(self.o74807, kmer_length=4)
        f.score(self.sprot12071)
        result = f.finish()
        self.assertEqual(len(result.indices[0]), 11899)
        self.assertEqual(
            result.indices[0][:10], 
            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
        )
        self.assertEqual(
            result.indices[0][-10:],
            [12061, 12062, 12063, 12064, 12065, 12066, 12067, 12068, 12069, 12070],
        )


class TestT50(_TestFilter, unittest.TestCase):
    @unittest.skipUnless(data.exists("uniprot_sprot15.fasta"), "missing data file")
    @unittest.skipUnless(data.exists("uniprot_sprot196.fasta"), "missing data file")
    def test_sprot15_sprot196(self):
        f = HeuristicFilter(self.sprot15, score_threshold=50)
        f.score(self.sprot196)
        result = f.finish()
        self.assertEqual(len(result.indices[0]), 76)
        self.assertEqual(len(result.indices[1]), 94)
        self.assertEqual(len(result.indices[2]), 107)
        self.assertEqual(len(result.indices[3]), 56)
        self.assertEqual(len(result.indices[4]), 113)
        self.assertEqual(len(result.indices[5]), 29)
        self.assertEqual(len(result.indices[6]), 62)
        self.assertEqual(len(result.indices[7]), 75)
        self.assertEqual(len(result.indices[8]), 111)
        self.assertEqual(len(result.indices[9]), 13)
        self.assertEqual(len(result.indices[10]), 54)
        self.assertEqual(len(result.indices[11]), 113)
        self.assertEqual(len(result.indices[12]), 79)
        self.assertEqual(len(result.indices[13]), 161)
        self.assertEqual(len(result.indices[14]), 24)

    @unittest.skipUnless(data.exists("O74807.fasta"), "missing data file")
    @unittest.skipUnless(data.exists("uniprot_sprot15.fasta"), "missing data file")
    def test_o74807_sprot15(self):
        f = HeuristicFilter(self.o74807, score_threshold=50)
        f.score(self.sprot15)
        result = f.finish()
        self.assertEqual(result.indices[0], [2, 8, 10, 12, 14])

    @unittest.skipUnless(data.exists("O74807.fasta"), "missing data file")
    @unittest.skipUnless(data.exists("uniprot_sprot196.fasta"), "missing data file")
    def test_o74807_sprot196(self):
        f = HeuristicFilter(self.o74807, score_threshold=50)
        f.score(self.sprot196)
        result = f.finish()
        self.assertEqual(len(result.indices[0]), 41)
        self.assertEqual(
            result.indices[0][:10],
            [2, 8, 10, 12, 14, 16, 18, 25, 26, 29],
        )
        self.assertEqual(
            result.indices[0][-10:],
            [156, 157, 158, 161, 164, 173, 174, 176, 180, 181]
        )

    @unittest.skipUnless(data.exists("O74807.fasta"), "missing data file")
    @unittest.skipUnless(data.exists("uniprot_sprot12071.fasta"), "missing data file")
    def test_o74807_sprot12071(self):
        f = HeuristicFilter(self.o74807, score_threshold=50)
        f.score(self.sprot12071)
        result = f.finish()
        self.assertEqual(len(result.indices[0]), 4444)
        self.assertEqual(
            result.indices[0][:10], 
            [2, 8, 10, 12, 14, 16, 18, 25, 26, 29],
        )
        self.assertEqual(
            result.indices[0][-10:],
            [12053, 12055, 12059, 12060, 12061, 12063, 12064, 12065, 12069, 12070],
        )