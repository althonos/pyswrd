# distutils: language = c++
# cython: language_level=3, linetrace=True, binding=True
"""Bindings to SWORD, a heuristic method for sequence database search.

References:
    - Korpar M., M. Šošić, D. Blažeka, M. Šikić.
      *SW#db: GPU-Accelerated Exact Sequence Similarity Database Search*.
      PLoS One. 2015;10(12):e0145857. Published 2015 Dec 31.
      :doi:`10.1371/journal.pone.0145857`.
    - Vaser R., D. Pavlović, M. Šikić. 
      *SWORD—a highly efficient protein database search*.
      Bioinformatics, Volume 32, Issue 17, September 2016, Pages i680–i684.
      :doi:`10.1093/bioinformatics/btw445`.

"""

from cython.operator cimport dereference, preincrement

cimport libcpp.algorithm
from libc.stdint cimport int32_t, uint16_t, uint32_t, uint64_t
from libcpp cimport bool, nullptr
from libcpp.pair cimport pair
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.memory cimport unique_ptr, shared_ptr

cimport sword.chain
cimport sword.hash
cimport sword.kmers
cimport sword.reader
cimport sword.score_matrix
from sword.database_search cimport ChainEntry as _ChainEntry, ChainEntrySet as _ChainEntrySet, Indexes as _Indexes
from sword.kmers cimport Kmers as _Kmers
from sword.chain cimport ChainSet as _ChainSet, Chain as _Chain
from sword.reader cimport Reader as _Reader 
from sword.hash cimport Iterator as _HashIterator
from sword.score_matrix cimport ScoreMatrix as _ScoreMatrix, ScoreMatrixType as _ScoreMatrixType

cimport pyopal.lib
from pyopal.lib cimport seq_t, digit_t

import os

# --- C++ constants and helpers -----------------------------------------------

cdef extern from * nogil:
    """
    bool chainLengthKey(const std::unique_ptr<Chain>& left, const std::unique_ptr<Chain>& right) {
        return left->length() < right->length();
    }
    """
    bool chainLengthKey(const unique_ptr[_Chain]& left, const unique_ptr[_Chain]& right) noexcept

cdef extern from * nogil:
    """
    bool chainEntryDataKey(const ChainEntry& left, const ChainEntry& right) {
        return left.data() > right.data();
    }
    """
    bool chainEntryDataKey(const _ChainEntry& left, const _ChainEntry& right) noexcept

cdef uint32_t    kProtBits   = 5
cdef uint32_t[6] kmerDelMask = [ 0, 0, 0, 0x7fff, 0xFFFFF, 0x1FFFFFF ]

# --- Python ------------------------------------------------------------------

cdef class Kmers:
    cdef shared_ptr[_Kmers] _kmers

    def __init__(self, ScoreMatrix score_matrix, kmer_length = 3, score_threshold = 13):
        self._kmers = shared_ptr[_Kmers](
            sword.kmers.createKmers(kmer_length, score_threshold, score_matrix._score_matrix)
        )


cdef class ChainSet:
    cdef _ChainSet _chains

    def __len__(self):
        return self._chains.size()

    def __getitem__(self, ssize_t i):
        cdef ssize_t _i = i
        if _i < 0:
            _i += self._chains.size()
        if _i < 0 or _i >= self._chains.size():
            raise IndexError(i)
        return self._chains[i].get().data()


cdef class Reader:
    cdef shared_ptr[_Reader] _reader

    def __init__(self, path):
        cdef bytes _path = os.fsencode(path)
        self._reader = shared_ptr[_Reader](sword.reader.createReader(_path))

    cpdef ChainSet read(self):
        assert self._reader != nullptr
        cdef ChainSet s = ChainSet()
        self._reader.get().read_chains(s._chains, 0)
        return s


cdef class ScoreMatrix:
    cdef shared_ptr[_ScoreMatrix] _score_matrix

    def __init__(self, str name = "BLOSUM62", int32_t gap_open = 10, int32_t gap_extend = 1):
        cdef _ScoreMatrixType ty
        if name == "BLOSUM62":
            ty = sword.score_matrix.kBlosum62
        else:
            raise ValueError(f"unsupported score matrix: {name!r}")
        self._score_matrix = shared_ptr[_ScoreMatrix](
            sword.score_matrix.createScoreMatrix(
                ty,
                gap_open,
                gap_extend,
            )
        )

    @property
    def name(self):
        assert self._score_matrix != nullptr
        return self._score_matrix.get().scorerName().decode()


cdef class FilterScore:
    cdef readonly uint32_t index
    cdef readonly uint32_t score 

    def __cinit__(self, uint32_t index, uint32_t score):
        self.index = index
        self.score = score

    def __repr__(self):
        return f"{type(self).__name__}({self.index!r}, {self.score!r})"


cdef class FilterResult:
    cdef readonly uint32_t database_size
    cdef readonly list     entries
    cdef readonly list     indices

    def __init__(self, database_size, entries, indices):
        self.entries = entries
        self.indices = indices
        self.database_size = database_size


cdef class HeuristicFilter:

    cdef readonly ChainSet       queries
    cdef readonly Kmers          kmers

    cdef readonly uint32_t       score_threshold
    cdef readonly uint32_t       max_candidates

    cdef          uint32_t       database_size
    cdef          _ChainEntrySet entries

    def __init__(self, ChainSet queries, *, kmer_length: uint32_t = 3, max_candidates: uint32_t = 30000, score_threshold = 13, ScoreMatrix score_matrix = ScoreMatrix()):
        self.queries = queries
        self.score_threshold = score_threshold
        self.max_candidates = max_candidates
        
        self.kmers = Kmers(score_matrix, kmer_length, score_threshold)
        
        self.database_size = 0
        self.entries = _ChainEntrySet(len(self.queries))

    cpdef void score(self, ChainSet database):
        cdef uint32_t      kmer
        cdef _HashIterator begin
        cdef _HashIterator end
        
        cdef uint32_t      database_size = database._chains.size()
        cdef uint32_t      database_start = 0
        cdef uint32_t      database_end   = database_size

        cdef _ChainEntrySet entries_part      = _ChainEntrySet(self.queries._chains.size())
        cdef vector[uint16_t] min_entry_score = vector[uint16_t](self.queries._chains.size())
        cdef vector[uint16_t] entries_found   = vector[uint16_t](self.queries._chains.size())

        cdef uint32_t id_
        for i in range(self.queries._chains.size()):
            id_ = self.queries._chains[i].get().id()
            entries_found[i] = self.entries[id_].size()
            min_entry_score[i] = 65000 if entries_found[i] == 0 else self.entries[id_].back().data()

        cdef uint32_t kmer_length       = self.kmers._kmers.get().kmer_length()
        cdef uint32_t max_target_length = dereference(libcpp.algorithm.max_element(  database._chains.begin() + database_start, database._chains.begin() + database_end, chainLengthKey )).get().length()
        cdef size_t   groups            = 0

        cdef uint32_t kmer_offset       = kmer_length - 1
        cdef uint32_t del_mask          = kmerDelMask[kmer_length]
        cdef uint32_t max_scores_length = 100000 if kmer_length == 3 else 500000

        cdef vector[uint16_t] scores = vector[uint16_t](max_scores_length)
        cdef vector[uint32_t] score_lengths = vector[uint32_t](self.queries._chains.size())
        cdef vector[uint32_t] score_starts = vector[uint32_t](self.queries._chains.size() + 1)
        cdef vector[uint16_t] max_score = vector[uint16_t](self.queries._chains.size())
        score_starts[0] = 0

        cdef uint32_t min_score = 1 if kmer_length == 3 else 0

        while i < self.queries._chains.size():
            groups += 1

            group_length = 0
            scores_length = 0

            for j in range(i, self.queries._chains.size()):

                length = self.queries._chains[j].get().length() + max_target_length - 2 * kmer_length + 1

                if scores_length + length > max_scores_length and group_length > 0:
                    break

                scores_length += length
                group_length += 1

            hash_ = sword.hash.createHash(self.queries._chains, i, group_length, self.kmers._kmers)

            for j in range(database_start, database_end):

                for k in range(group_length):
                    score_lengths[k] = self.queries._chains[i + k].get().length() + database._chains[j].get().length() - 2 * kmer_length + 1
                    score_starts[k + 1] = score_starts[k] + score_lengths[k]

                sequence = database._chains[j].get().data()
                kmer = sequence[0]

                for k in range(1, kmer_offset):
                    kmer = (kmer << kProtBits) | sequence[k]

                max_diag_id = database._chains[j].get().length() - kmer_length
                for k in range(kmer_offset, sequence.size()):
                    kmer = ((kmer << kProtBits) | sequence[k]) & del_mask
                    hash_.get().hits(begin, end, kmer)

                    while begin != end:
                        diagonal = max_diag_id + kmer_offset - k + dereference(begin).position() + score_starts[dereference(begin).id()]
                        scores[diagonal] += 1
                        if max_score[dereference(begin).id()] < scores[diagonal]:
                            max_score[dereference(begin).id()] = scores[diagonal]
                        preincrement(begin)

                for k in range(group_length):
                    if max_score[k] <= min_score:
                        continue

                    id_ = self.queries._chains[i + k].get().id()
                    flag = entries_part[id_].size() < self.max_candidates and entries_found[k] < self.max_candidates

                    if flag or max_score[k] >= min_entry_score[k]:
                        entries_part[id_].emplace_back(_ChainEntry(database._chains[j].get().id(), max_score[k]))
                        if min_entry_score[k] > max_score[k]:
                            min_entry_score[k] = max_score[k]

                for k in range(group_length):
                    if max_score[k] == 0:
                        continue
                    
                    max_score[k] = 0
                    libcpp.algorithm.fill_n(&scores[0] + score_starts[k], score_lengths[k], 0)

            for k in range(group_length):
                id_ = self.queries._chains[i + k].get().id()
                self.entries[id_].insert( self.entries[id_].end(), entries_part[id_].begin(), entries_part[id_].end() )
                entries_part[id_].clear()

                libcpp.algorithm.stable_sort(self.entries[id_].begin(), self.entries[id_].end(), chainEntryDataKey);

                if self.entries[id_].size() > self.max_candidates:
                    self.entries[id_].resize(self.max_candidates)

            i += group_length
            print(i)

        self.database_size += database_size

    cpdef FilterResult finish(self):
        entries = [
            [FilterScore(entry.chain_idx(), entry.data()) for entry in entries]
            for entries in self.entries
        ]
        indices = [
            sorted([entry.index for entry in x]) for x in entries
        ]
        return FilterResult(entries=entries, indices=indices, database_size=self.database_size)



