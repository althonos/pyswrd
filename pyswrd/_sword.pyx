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
cimport libcpp.utility
from libc.stdint cimport int32_t, uint16_t, uint32_t, uint64_t
from libcpp cimport bool, nullptr
from libcpp.pair cimport pair
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp.memory cimport unique_ptr, shared_ptr
from mutex cimport mutex, unique_lock

cimport sword.chain
cimport sword.hash
cimport sword.kmers
cimport sword.reader
cimport sword.score_matrix
from sword.database_search cimport ChainEntry as _ChainEntry, ChainEntrySet as _ChainEntrySet, Indexes as _Indexes
from sword.kmers cimport Kmers as _Kmers
from sword.chain cimport ChainSet as _ChainSet, Chain as _Chain
from sword.hash cimport Iterator as _HashIterator
from sword.score_matrix cimport ScoreMatrix as _ScoreMatrix, ScoreMatrixType as _ScoreMatrixType

import os
import functools
import multiprocessing.pool

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

    template<typename T>
    void stable_sort_by_length(T begin, T end) {
        std::stable_sort(begin, end, chainEntryDataKey);
    }
    """
    bool chainEntryDataKey(const _ChainEntry& left, const _ChainEntry& right) noexcept
    void stable_sort_by_length[T](T begin, T end)

cdef uint32_t    kProtBits   = 5
cdef uint32_t[6] kmerDelMask = [ 0, 0, 0, 0x7fff, 0xFFFFF, 0x1FFFFFF ]

# --- Python ------------------------------------------------------------------

cdef class KmerGenerator:
    """A generator of k-mers with optional substitutions.
    """
    cdef shared_ptr[_Kmers] _kmers

    def __init__(self, Scorer scorer, size_t kmer_length = 3, size_t score_threshold = 13):
        if kmer_length < 3 or kmer_length > 5:
            raise ValueError(f"kmer_length must be 3, 4 or 5, got: {kmer_length!r}")
        self._kmers = shared_ptr[_Kmers](
            sword.kmers.createKmers(kmer_length, score_threshold, scorer._sm)
        )


cdef class Sequences:
    """A list of sequences.
    """
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

    def append(self, str sequence):
        cdef bytes    seq   = sequence.encode('ascii')
        cdef uint32_t total = self._chains.size()
        cdef bytes    name  = str(total).encode()
        cdef unique_ptr[_Chain] chain = sword.chain.createChain( total, name, len(name), seq, len(seq) )
        self._chains.push_back(libcpp.utility.move(chain))


cdef class Scorer:
    """A class storing the scoring matrix and gap parameters for alignments.
    """
    cdef shared_ptr[_ScoreMatrix] _sm

    def __init__(self, str name = "BLOSUM62", int32_t gap_open = 10, int32_t gap_extend = 1):
        cdef _ScoreMatrixType ty
        if name == "BLOSUM62":
            ty = sword.score_matrix.kBlosum62
        else:
            raise ValueError(f"unsupported score matrix: {name!r}")
        self._sm = shared_ptr[_ScoreMatrix](
            sword.score_matrix.createScoreMatrix(
                ty,
                gap_open,
                gap_extend,
            )
        )

    @property
    def gap_open(self):
        """`int`: The score penalty for creating a gap.
        """
        return self._sm.gap_open()

    @property
    def gap_extend(self):
        """`int`: The score penalty for extending a gap.
        """
        return self._sm.gap_extend()

    @property
    def name(self):
        """`str`: The name of the scoring matrix.
        """
        assert self._sm != nullptr
        return self._sm.get().scorerName().decode('ascii')


cdef class FilterScore:
    """The score of the heuristic filter for a single target. 

    Attributes:
        index (`int`): The index of the sequence in the target database.
        score (`int`): The score of the sequence.

    """
    cdef readonly uint32_t index
    cdef readonly uint32_t score 

    def __cinit__(self, uint32_t index, uint32_t score):
        self.index = index
        self.score = score

    def __repr__(self):
        return f"{type(self).__name__}({self.index!r}, {self.score!r})"


cdef class FilterResult:
    """The result of the heuristic filter.
    """
    cdef readonly uint32_t database_size
    cdef readonly list     entries
    cdef readonly list     indices

    def __init__(self, uint32_t database_size, list entries, list indices):
        self.entries = entries
        self.indices = indices
        self.database_size = database_size


cdef class HeuristicFilter:
    """The SWORD heuristic filter for selecting alignment candidates.
    """
    cdef readonly Sequences                 queries
    cdef readonly KmerGenerator             kmers
    cdef readonly uint32_t                  score_threshold

    cdef readonly uint32_t                  max_candidates
    cdef          uint32_t                  database_size
    cdef          _ChainEntrySet            entries

    cdef readonly size_t                    threads
    cdef readonly object                    pool
    cdef          vector[unique_ptr[mutex]] mutexes

    def __init__(self, Sequences queries, *, kmer_length: uint32_t = 3, max_candidates: uint32_t = 30000, score_threshold = 13, Scorer scorer = Scorer(), size_t threads = 0):
        # k-mer generation parameter
        self.queries = queries
        self.score_threshold = score_threshold
        self.kmers = KmerGenerator(scorer, kmer_length, score_threshold)
        # parameters and buffers for candidate retrieval
        self.max_candidates = max_candidates
        self.database_size = 0
        self.entries = _ChainEntrySet(len(self.queries))
        # parameters for multiprocessing
        self.threads = threads or os.cpu_count()
        self.pool = multiprocessing.pool.ThreadPool(self.threads)
        self.mutexes = vector[unique_ptr[mutex]]()
        for i in range(len(self.queries)):
            self.mutexes.push_back( unique_ptr[mutex]( new mutex() ) )

    def __del__(self):
        self.pool.close()
        self.pool.join()

    cpdef vector[uint32_t] _preprocess_database(
        self,
        Sequences database,
        uint32_t max_short_length = 2000,
    ):
        # NOTE: ported from: preprocDatabase in `database_search.cpp`
        # FIXME: at the moment doesn't work properly because the chains are 
        #        not sorted by length.

        cdef vector[uint32_t] dst
        cdef uint32_t         i
        cdef uint64_t         total_length = 0
        cdef uint64_t         short_total_length = 0
        cdef uint64_t         long_total_length  = 0
        cdef uint32_t         split              = 0
        cdef uint64_t         short_task_size
        cdef uint64_t         long_task_size

        # sort by length (FIXME?)
        # libcpp.algorithm.sort(database._chains.begin(), database._chains.end(), chainLengthKey)

        # split tasks between long and short
        for i in range(database._chains.size()):
            l = database._chains[i].get().length() 
            if l > max_short_length:
                if split == 0:
                    split = i
                long_total_length += l
            else:
                short_total_length += l

        if short_total_length == 0:
            split = 0
        if long_total_length == 0:
            split = database._chains.size()

        # spread tasks across threads
        short_task_size = short_total_length / self.threads
        long_task_size = long_total_length / self.threads

        dst.reserve(2*self.threads + 1)
        dst.emplace_back(0)

        total_length = 0
        for i in range(split):
            total_length += database._chains[i].get().length()
            if total_length > short_task_size:
                total_length = 0
                dst.emplace_back(i + 1)
        
        if dst.back() != split:
            dst.emplace_back(split)

        total_length = 0
        for i in range(split, database._chains.size()):
            total_length += database._chains[i].get().length()
            if total_length > long_task_size:
                total_length = 0
                dst.emplace_back(i + 1)

        if dst.back() != database._chains.size():
            dst.emplace_back(database._chains.size())

        return dst

    cpdef void _score_chunk(self, Sequences database, uint32_t database_start, uint32_t database_end):
        cdef uint32_t      kmer
        cdef _HashIterator begin
        cdef _HashIterator end
        cdef uint32_t      i
        cdef uint32_t      j
        cdef uint32_t      k
        cdef uint32_t      diagonal
        cdef uint32_t      max_diag_id
        cdef bool          flag

        cdef uint32_t      database_size = database._chains.size()

        cdef unique_lock[mutex] lock

        cdef _ChainEntrySet entries_part      = _ChainEntrySet(self.queries._chains.size())
        cdef vector[uint16_t] min_entry_score = vector[uint16_t](self.queries._chains.size())
        cdef vector[uint16_t] entries_found   = vector[uint16_t](self.queries._chains.size())

        cdef uint32_t id_
        with nogil:
            for i in range(self.queries._chains.size()):
                id_ = self.queries._chains[i].get().id()
                lock = unique_lock[mutex](self.mutexes[id_].get()[0])
                entries_found[i] = self.entries[id_].size()
                min_entry_score[i] = 65000 if entries_found[i] == 0 else self.entries[id_].back().data()
                lock.unlock()

        cdef uint32_t length
        cdef uint32_t kmer_length       = self.kmers._kmers.get().kmer_length()
        cdef uint32_t max_target_length = dereference(libcpp.algorithm.max_element( database._chains.begin() + database_start, database._chains.begin() + database_end, chainLengthKey )).get().length()
        cdef size_t   groups            = 0
        cdef size_t   group_length      = 0

        cdef uint32_t kmer_offset       = kmer_length - 1
        cdef uint32_t del_mask          = kmerDelMask[kmer_length]
        cdef uint32_t scores_length     = 0
        cdef uint32_t max_scores_length = 100000 if kmer_length == 3 else 500000

        cdef vector[uint16_t] scores = vector[uint16_t](max_scores_length)
        cdef vector[uint32_t] score_lengths = vector[uint32_t](self.queries._chains.size())
        cdef vector[uint32_t] score_starts = vector[uint32_t](self.queries._chains.size() + 1)
        cdef vector[uint16_t] max_score = vector[uint16_t](self.queries._chains.size())
        score_starts[0] = 0

        cdef uint32_t min_score = 1 if kmer_length == 3 else 0

        with nogil:

            i = 0
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

                    max_diag_id = database._chains[j].get().length() - kmer_length
                    sequence = database._chains[j].get().data()
                    kmer = sequence[0]

                    for k in range(1, kmer_offset):
                        kmer = (kmer << kProtBits) | sequence[k]
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
                    lock = unique_lock[mutex](self.mutexes[id_].get()[0])
                    self.entries[id_].insert( self.entries[id_].end(), entries_part[id_].begin(), entries_part[id_].end() )
                    entries_part[id_].clear()
                    stable_sort_by_length(self.entries[id_].begin(), self.entries[id_].end())
                    if self.entries[id_].size() > self.max_candidates:
                        self.entries[id_].resize(self.max_candidates)
                    lock.unlock()

                i += group_length

    cpdef void score(self, Sequences database):
        if self.threads > 1:
            splits = list(self._preprocess_database(database))
            score = functools.partial(self._score_chunk, database)
            self.pool.starmap(score,  zip(splits[:-1], splits[1:]))
        else:
            self._score_chunk(database, 0, len(database))
        self.database_size += len(database)
        
    cpdef FilterResult finish(self):

        self.pool.close()

        entries = [
            [FilterScore(entry.chain_idx(), entry.data()) for entry in entries]
            for entries in self.entries
        ]
        indices = [
            sorted([entry.index for entry in x]) for x in entries
        ]
        return FilterResult(entries=entries, indices=indices, database_size=self.database_size)


    
