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
from cpython.buffer cimport PyBUF_FORMAT, PyBUF_READ, PyBUF_WRITE

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
cimport sword.evalue
cimport sword.score_matrix
from sword.database_search cimport ChainEntry as _ChainEntry, ChainEntrySet as _ChainEntrySet, Indexes as _Indexes
from sword.kmers cimport Kmers as _Kmers
from sword.chain cimport ChainSet as _ChainSet, Chain as _Chain
from sword.hash cimport Iterator as _HashIterator
from sword.evalue cimport EValue as _EValue
from sword.score_matrix cimport ScoreMatrix as _ScoreMatrix, ScoreMatrixType as _ScoreMatrixType

cimport pyopal.lib
from pyopal.lib cimport digit_t

import contextlib
import operator
import os
import functools
import threading
from string import ascii_uppercase
from pyopal import align

try:
    from multiprocessing.pool import ThreadPool
except ImportError:
    ThreadPool = None  # multiprocessing.pool may be missing on some environments


# --- C++ constants and helpers -----------------------------------------------

cdef extern from * nogil:
    """
    template<typename T1, typename T2>
    void stable_sort_by_second(std::vector<std::pair<T1, T2>>& v) {
        std::stable_sort(v.begin(), v.end(), [](const auto& a, const auto& b) {return a.second < b.second; });
        // std::stable_sort(begin, end, pairSecondKey);
    }
    """
    void stable_sort_by_second[T1, T2](vector[pair[T1, T2]]& v)

cdef extern from * nogil:
    """
    bool chainLengthKey(const std::shared_ptr<Chain>& left, const std::shared_ptr<Chain>& right) {
        return left->length() < right->length();
    }
    """
    bool chainLengthKey(const shared_ptr[_Chain]& left, const shared_ptr[_Chain]& right) noexcept

cdef extern from * nogil:
    """
    void stable_sort_by_length(std::vector<ChainEntry>& v) {
        std::stable_sort(v.begin(), v.end(), [](const auto& a, const auto& b) {return (a.data(), a.chain_idx()) > (b.data(), b.chain_idx()); });
    }
    """
    void stable_sort_by_length(vector[_ChainEntry] v)

cdef uint32_t    kProtBits   = 5
cdef uint32_t[6] kmerDelMask = [ 0, 0, 0, 0x7fff, 0xFFFFF, 0x1FFFFFF ]

# --- Constants ----------------------------------------------------------------

cdef Py_ssize_t[2] _SWORD_SCORE_MATRIX_SHAPE = [
    sword.score_matrix.num_rows_,
    sword.score_matrix.num_columns_,
]

cdef dict _SWORD_SCORE_MATRICES = {
    "BLOSUM45": <int> sword.score_matrix.kBlosum45,
    "BLOSUM50": <int> sword.score_matrix.kBlosum50,
    "BLOSUM62": <int> sword.score_matrix.kBlosum62,
    "BLOSUM80": <int> sword.score_matrix.kBlosum80,
    "BLOSUM90": <int> sword.score_matrix.kBlosum90,
    "PAM30": <int> sword.score_matrix.kPam30,
    "PAM70": <int> sword.score_matrix.kPam70,
    "PAM250": <int> sword.score_matrix.kPam250,
}

cdef pyopal.lib.Alphabet _SWORD_ALPHABET = pyopal.lib.Alphabet(
    ascii_uppercase
)

# --- Python helpers -----------------------------------------------------------

@contextlib.contextmanager
def nullcontext(enter_result):
    """Return a context manager that returns its input and does nothing.

    Adapted from `contextlib.nullcontext` for backwards compatibility
    with Python 3.6.

    """
    yield enter_result

# --- Parameters ---------------------------------------------------------------

cdef class KmerGenerator:
    """A generator of k-mers with optional substitutions.
    """
    cdef          shared_ptr[_Kmers] _kmers
    cdef readonly Scorer             scorer

    def __init__(self, Scorer scorer, size_t kmer_length = 3, size_t score_threshold = 13):
        self.scorer = scorer
        if kmer_length < 3 or kmer_length > 5:
            raise ValueError(f"kmer_length must be 3, 4 or 5, got: {kmer_length!r}")
        self._kmers = shared_ptr[_Kmers](
            sword.kmers.createKmers(kmer_length, score_threshold, scorer._sm)
        )

cdef class Scorer:
    """A class storing the scoring matrix and gap parameters for alignments.
    """
    cdef          shared_ptr[_ScoreMatrix] _sm
    cdef readonly pyopal.lib.ScoreMatrix   score_matrix

    def __init__(self, str name = "BLOSUM62", int32_t gap_open = 10, int32_t gap_extend = 1):
        cdef _ScoreMatrixType ty
        if name in _SWORD_SCORE_MATRICES:
            ty = <_ScoreMatrixType> <int> _SWORD_SCORE_MATRICES[name]
        else:
            raise ValueError(f"unsupported score matrix: {name!r}")
        self._sm = shared_ptr[_ScoreMatrix](
            sword.score_matrix.createScoreMatrix(
                ty,
                gap_open,
                gap_extend,
            )
        )

        cdef int* scores = self._sm.get().data()
        self.score_matrix = pyopal.lib.ScoreMatrix(
            alphabet=_SWORD_ALPHABET,
            matrix = [
                [
                    scores[i*sword.score_matrix.num_columns_ + j]
                    for j in range( sword.score_matrix.num_columns_ )
                ]
                for i in range(sword.score_matrix.num_rows_)
            ]
        )

    def __repr__(self):
        cdef str ty = type(self).__name__
        return f"{ty}({self.name!r}, gap_open={self.gap_open!r}, gap_extend={self.gap_extend!r})"

    def __reduce__(self):
        return (type(self), (self.name, self.gap_open, self.gap_extend))

    def __getbuffer__(self, Py_buffer* buffer, int flags):
        if flags & PyBUF_FORMAT:
            buffer.format = b"i"
        else:
            buffer.format = NULL
        buffer.buf = self._sm.get().data()
        buffer.internal = NULL
        buffer.itemsize = sizeof(int)
        buffer.len = sword.score_matrix.num_rows_ * sword.score_matrix.num_columns_ * sizeof(int)
        buffer.ndim = 2
        buffer.obj = self
        buffer.readonly = 1
        buffer.shape = <Py_ssize_t*> _SWORD_SCORE_MATRIX_SHAPE
        buffer.suboffsets = NULL
        buffer.strides = NULL

    @property
    def gap_open(self):
        """`int`: The score penalty for creating a gap.
        """
        assert self._sm != nullptr
        return self._sm.get().gap_open()

    @property
    def gap_extend(self):
        """`int`: The score penalty for extending a gap.
        """
        assert self._sm != nullptr
        return self._sm.get().gap_extend()

    @property
    def name(self):
        """`str`: The name of the scoring matrix.
        """
        assert self._sm != nullptr
        return self._sm.get().scorerName().decode('ascii')

cdef class EValue:
    """A class for calculating E-values from alignment scores.
    """
    cdef readonly Scorer              scorer
    cdef          shared_ptr[_EValue] _evalue

    def __init__(self, uint64_t database_size, Scorer scorer):
        self.scorer = scorer
        self._evalue = shared_ptr[_EValue](sword.evalue.createEValue(database_size, scorer._sm))

    cpdef double calculate(self, int32_t score, uint32_t query_length, uint32_t target_length):
        return self._evalue.get().calculate(score, query_length, target_length)

# --- Sequence Storage ---------------------------------------------------------

ctypedef fused Indices:
    object
    vector[uint32_t]

ctypedef fused Mask:
    object
    vector[bool]

cdef class Sequences(pyopal.lib.BaseDatabase):
    """A list of sequences.
    """
    cdef _ChainSet        _chains
    cdef vector[digit_t*] _pointers
    cdef vector[int]      _lengths

    # --- Magic methods --------------------------------------------------------

    def __cinit__(self):
        self._chains = _ChainSet()

    def __init__(self, object sequences=()):
        super().__init__(alphabet=_SWORD_ALPHABET)
        self.clear()
        self.extend(sequences)

    def __reduce__(self):
        return (type(self), ((),), None, iter(self))

    # --- Database interface ---------------------------------------------------

    cdef size_t get_size(self) noexcept:
        return self._chains.size()

    cdef const digit_t** get_sequences(self) except NULL:
        cdef digit_t** sequences = self._pointers.data()
        return <const digit_t**> sequences

    cdef const int* get_lengths(self) except NULL:
        return self._lengths.data()

    # --- Sequence interface ---------------------------------------------------

    def __getitem__(self, object index):
        if isinstance(index, slice):
            with self.lock.read:
                size = self.get_size()
                indices = range(*index.indices(size))
                return self.extract(indices)
        return super().__getitem__(index)

    cpdef void clear(self) except *:
        """Remove all sequences from the database.
        """
        with self.lock.write:
            self._chains.clear()
            self._pointers.clear()
            self._lengths.clear()

    cpdef void append(self, object sequence):
        """Append a single sequence at the end of the sequence list.

        Arguments:
            sequence (`str` or byte-like object): The new sequence.

        Hint:
            When inserting several sequences in the list, consider
            using the `Sequences.extend` method instead so that the
            internal buffers can reserve space just once for every
            new sequence.

        Example:
            >>> db = pyswrd.Sequences(["ATGC", "TTCA"])
            >>> db.append("AAAA")
            >>> list(db)
            ['ATGC', 'TTCA', 'AAAA']

        """
        cdef bytes              seq
        cdef shared_ptr[_Chain] chain
        cdef uint32_t           total = self._chains.size()
        cdef bytes              name  = str(total).encode()

        seq = sequence.encode('ascii') if isinstance(sequence, str) else sequence
        chain = shared_ptr[_Chain](sword.chain.createChain(total, name, len(name), seq, len(seq)))

        with self.lock.write:
            self._lengths.push_back(chain.get().length())
            self._pointers.push_back(<digit_t*> chain.get().data().data())
            self._chains.push_back(chain)

    cpdef void extend(self, object sequences) except *:
        """Extend the sequence list by adding sequences from an iterable.

        Example:
            >>> db = pyswrd.Sequences(["ATGC"])
            >>> db.extend(["TTCA", "AAAA", "GGTG"])
            >>> list(db)
            ['ATGC', 'TTCA', 'AAAA', 'GGTG']

        """
        cdef size_t size
        cdef size_t hint = operator.length_hint(sequences)

        # attempt to reserve space in advance for the new sequences
        with self.lock.write:
            size = self._chains.size()
            if hint > 0:
                self._chains.reserve(hint + size)
                self._pointers.reserve(hint + size)
                self._lengths.reserve(hint + size)

        # append sequences in order
        for sequence in sequences:
            self.append(sequence)

    cpdef void reverse(self) except *:
        raise NotImplementedError("Sequences.reverse")

    cpdef void insert(self, ssize_t index, object sequence):
        raise NotImplementedError("Sequences.insert")

    cpdef Sequences mask(self, Mask bitmask):
        cdef bool      b
        cdef size_t    i
        cdef Sequences sequences
        cdef size_t    mask_size

        sequences = Sequences.__new__(Sequences)
        sequences.alphabet = self.alphabet

        with self.lock.read:
            if Mask is object:
                mask_size = len(bitmask)
            else:
                mask_size = bitmask.size()

            if mask_size != len(self):
                raise IndexError(bitmask)

            with nogil(Mask is not object):
                i = 0
                for b in bitmask:
                    if b:
                        sequences._chains.push_back(self._chains[i])
                        sequences._pointers.push_back(self._pointers[i])
                        sequences._lengths.push_back(self._lengths[i])
                    i += 1

        return sequences

    cpdef Sequences extract(self, Indices indices):
        cdef size_t    i
        cdef size_t    indices_size
        cdef Sequences sequences
        cdef size_t    size

        sequences = Sequences.__new__(Sequences)
        sequences.alphabet = self.alphabet

        with self.lock.read:
            # get length of database and indices
            size = self.get_size()
            if Indices is object:
                indices_size = len(indices)
            else:
                indices_size = indices.size()
            # reserve exact amount of memory
            with nogil:
                sequences._chains.reserve(indices_size)
                sequences._pointers.reserve(indices_size)
                sequences._lengths.reserve(indices_size)
            # recover sequences by index
            if Indices is object:
                for i in indices:
                    if i < 0 or i >= size:
                        raise IndexError(i)
                    sequences._chains.push_back(self._chains[i])
                    sequences._pointers.push_back(self._pointers[i])
                    sequences._lengths.push_back(self._lengths[i])
            else:
                with nogil:
                    for i in indices:
                        if i < 0 or i >= size:
                            raise IndexError(i)
                        sequences._chains.push_back(self._chains[i])
                        sequences._pointers.push_back(self._pointers[i])
                        sequences._lengths.push_back(self._lengths[i])

        return sequences

# --- Heuristic Filter ---------------------------------------------------------

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

    def __eq__(self, object other):
        if not isinstance(other, FilterScore):
            return NotImplemented
        return self.index == other.index and self.score == other.score


cdef class FilterResult:
    """The result of the heuristic filter.
    """
    cdef readonly uint32_t database_size
    cdef readonly uint64_t database_length
    cdef readonly list     entries
    cdef readonly list     indices
    cdef          _Indexes _indices

    def __init__(self, uint32_t database_size, uint64_t database_length, list entries, list indices):
        self.entries = entries
        self._indices = self.indices = indices
        self.database_size = database_size
        self.database_length = database_length


cdef class HeuristicFilter:
    """The SWORD heuristic filter for selecting alignment candidates.
    """
    cdef readonly Sequences                 queries
    cdef readonly KmerGenerator             kmer_generator
    cdef readonly uint32_t                  score_threshold

    cdef readonly uint32_t                  max_candidates
    cdef          uint32_t                  database_size
    cdef          uint64_t                  database_length
    cdef          _ChainEntrySet            entries

    cdef readonly size_t                    threads
    cdef readonly object                    pool
    cdef          vector[unique_ptr[mutex]] mutexes

    def __init__(
        self,
        Sequences queries not None,
        *,
        kmer_length: uint32_t = 3,
        max_candidates: uint32_t = 30000,
        score_threshold = 13,
        Scorer scorer = Scorer(),
        size_t threads = 0,
    ):
        # k-mer generation parameter
        self.queries = queries
        self.score_threshold = score_threshold
        self.kmer_generator = KmerGenerator(scorer, kmer_length, score_threshold)
        # parameters and buffers for candidate retrieval
        self.max_candidates = max_candidates
        self.database_size = 0
        self.database_length = 0
        self.entries = _ChainEntrySet(len(self.queries))
        # parameters for multiprocessing
        self.threads = os.cpu_count() if threads == 0 else threads
        if ThreadPool is None:
            self.threads = 1
        elif self.threads > 1:
            self.pool = ThreadPool(self.threads)
        self.mutexes = vector[unique_ptr[mutex]]()
        for i in range(len(self.queries)):
            self.mutexes.push_back(unique_ptr[mutex](new mutex()))

    @property
    def scorer(self):
        """`pyswrd.Scorer`: The scorer used for generating the k-mers.
        """
        return self.kmer_generator.scorer

    cpdef vector[uint32_t] _preprocess_database_long_short(
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

        # sort by length (FIXME!)
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

    cpdef vector[uint32_t] _preprocess_database(
        self,
        Sequences database,
    ):
        cdef vector[uint32_t] dst
        cdef uint32_t         i
        cdef uint32_t         task_size
        cdef uint32_t         task_length  = 0
        cdef uint64_t         total_length = 0
        cdef uint32_t         split        = 0

        # compute total length of database sequences
        for i in range(database._chains.size()):
            l = database._chains[i].get().length()
            total_length += l

        # spread tasks equally across threads
        task_size = total_length / self.threads
        dst.reserve(self.threads + 1)
        dst.emplace_back(0)

        # build groups of equal size
        for i in range(database._chains.size()):
            task_length += database._chains[i].get().length()
            if task_length > task_size:
                task_length = 0
                dst.emplace_back(i + 1)

        # make sure all the database is covered
        if dst.back() != database._chains.size():
            dst.emplace_back(database._chains.size())

        return dst

    cpdef void _score_chunk(self, Sequences database, uint32_t database_start, uint32_t database_end):
        cdef uint32_t           kmer
        cdef _HashIterator      begin
        cdef _HashIterator      end
        cdef uint32_t           i
        cdef uint32_t           j
        cdef uint32_t           k
        cdef uint32_t           id_
        cdef uint32_t           diagonal
        cdef uint32_t           max_diag_id
        cdef uint32_t           database_id
        cdef bool               flag
        cdef unique_lock[mutex] lock

        cdef vector[uint16_t]   scores
        cdef vector[uint32_t]   score_lengths
        cdef vector[uint32_t]   score_starts
        cdef vector[uint16_t]   max_score

        cdef uint32_t           length
        cdef uint32_t           max_target_length
        cdef uint32_t           queries_size      = self.queries._chains.size()
        cdef uint64_t           database_size     = database._chains.size()
        cdef _ChainEntrySet     entries_part      = _ChainEntrySet(queries_size)
        cdef vector[uint16_t]   min_entry_score   = vector[uint16_t](queries_size)
        cdef vector[uint16_t]   entries_found     = vector[uint16_t](queries_size)
        cdef uint32_t           kmer_length       = self.kmer_generator._kmers.get().kmer_length()
        cdef size_t             groups            = 0
        cdef size_t             group_length      = 0
        cdef uint32_t           kmer_offset       = kmer_length - 1
        cdef uint32_t           del_mask          = kmerDelMask[kmer_length]
        cdef uint32_t           scores_length     = 0
        cdef uint32_t           max_scores_length = 100000 if kmer_length == 3 else 500000
        cdef uint32_t           min_score         = 1 if kmer_length == 3 else 0

        with nogil:
            # Record the minimum score required for each query
            for i in range(queries_size):
                id_ = self.queries._chains[i].get().id()
                lock = unique_lock[mutex](self.mutexes[id_].get()[0])
                entries_found[i] = self.entries[id_].size()
                min_entry_score[i] = 65000 if entries_found[i] == 0 else self.entries[id_].back().data()
                lock.unlock()

            # Allocate space to store the scores
            scores = vector[uint16_t](max_scores_length)
            score_lengths = vector[uint32_t](queries_size)
            score_starts = vector[uint32_t](queries_size + 1)
            max_score = vector[uint16_t](queries_size)
            score_starts[0] = 0

            # Find the largest target sequence in the database chunk
            max_target_length = (
                dereference(
                    libcpp.algorithm.max_element(
                        database._chains.begin() + database_start,
                        database._chains.begin() + database_end,
                        chainLengthKey
                    )
                )
                .get()
                .length()
            )

            # Process all queries
            i = 0
            while i < queries_size:

                groups += 1
                group_length = 0
                scores_length = 0

                # Compute how many queries can be processed at the same time
                for j in range(i, queries_size):
                    length = self.queries._chains[j].get().length() + max_target_length - 2 * kmer_length + 1
                    if scores_length + length > max_scores_length and group_length > 0:
                        break
                    scores_length += length
                    group_length += 1

                # Compute hashes for the query group
                hash_ = sword.hash.createHash(self.queries._chains, i, group_length, self.kmer_generator._kmers)

                # Compare targets to the hashes
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
                            database_id = j + self.database_size
                            entries_part[id_].emplace_back(_ChainEntry(database_id, max_score[k]))
                            if min_entry_score[k] > max_score[k]:
                                min_entry_score[k] = max_score[k]

                    for k in range(group_length):
                        if max_score[k] == 0:
                            continue
                        max_score[k] = 0
                        libcpp.algorithm.fill_n(&scores[0] + score_starts[k], score_lengths[k], 0)

                # Merge the entries found in the query group with the entries found previously
                for k in range(group_length):
                    id_ = self.queries._chains[i + k].get().id()
                    lock = unique_lock[mutex](self.mutexes[id_].get()[0])
                    self.entries[id_].insert( self.entries[id_].end(), entries_part[id_].begin(), entries_part[id_].end() )
                    entries_part[id_].clear()
                    stable_sort_by_length(self.entries[id_])
                    if self.entries[id_].size() > self.max_candidates:
                        self.entries[id_].resize(self.max_candidates)
                    lock.unlock()

                # Advance to the next group
                i += group_length

    cpdef HeuristicFilter score(self, Sequences database):
        if self.threads > 1:
            splits = list(self._preprocess_database(database))
            score = functools.partial(self._score_chunk, database)
            self.pool.starmap(score,  zip(splits, splits[1:]))
        else:
            self._score_chunk(database, 0, len(database))
        self.database_size += len(database)
        for l in database._lengths:
            self.database_length += l
        return self

    cpdef FilterResult finish(self):
        if self.pool is not None:
            self.pool.close()
        entries = [
            [FilterScore(entry.chain_idx(), entry.data()) for entry in entries]
            for entries in self.entries
        ]
        indices = [
            sorted([entry.index for entry in x]) for x in entries
        ]
        return FilterResult(entries=entries, indices=indices, database_size=self.database_size, database_length=self.database_length)

# --- Database Search ---------------------------------------------------------

cdef class Hit:
    cdef readonly uint32_t              query_index
    cdef readonly uint32_t              target_index
    cdef readonly double                evalue
    cdef readonly pyopal.lib.FullResult result

    def __init__(self, query_index, target_index, evalue, result):
        self.query_index = query_index
        self.target_index = target_index
        self.evalue = evalue
        self.result = result

    def __repr__(self):
        cdef str ty = type(self).__name__
        return f"{ty}(query_index={self.query_index!r}, target_index={self.target_index!r}, evalue={self.evalue!r}, result={self.result!r})"

    @property
    def score(self):
        return self.result.score

def search(
    queries,
    targets,
    *,
    uint32_t gap_open = 10,
    uint32_t gap_extend = 1,
    str scorer_name = "BLOSUM62",
    uint32_t kmer_length = 3,
    uint32_t max_candidates = 30000,
    uint32_t score_threshold = 13,
    uint32_t max_alignments = 10,
    double max_evalue = 10.0,
    str algorithm = "sw",
    uint32_t threads = 0,
):
    """Run a many-to-many search of query sequences to target sequences.

    This function is a high-level wrapper around the different classes of
    the `pyswrd` library to support fast searches when all sequences are
    in memory.

    Arguments:
        queries (`~pyswrd.Sequences`, or iterable of `str`): The sequences
            to query the target sequences with.
        targets (`~pyswrd.Sequences` or iterable of `str`): The sequences
            to be queries with the query sequences.
        gap_open (`int`): The penalty for opening a gap in each alignment.
        gap_extend (`int`): The penalty for extending a gap in each
            alignment.
        scorer_name (`str`): The name of the scoring matrix to use
            for scoring each alignment. See `~pyswrd.Scorer` for the
            list of supported names.
        kmer_length (`int`): The length of the k-mers to use in the
            SWORD heuristic filter.
        max_candidates (`int`): The maximum number of candidates to
            retain in the heuristic filter.
        max_evalue (`float`): The E-value threshold above which to
            discard sequences before alignment.
        algorithm (`str`): The algorithm to use to perform pairwise
            alignment. See `pyopal.Database.search` for more
            information.
        threads (`int`): The number of threads to use to run the
            pre-filter and alignments. If zero is given, uses the
            number of CPUs reported by `os.cpu_count`. If one given,
            use the main threads for aligning, otherwise spawns a
            `multiprocessing.pool.ThreadPool`.

    Yields:
        `~pyswrd.Hit`: Hit objects for each hit passing the threshold
        parameters. Hits are grouped by query index, and sorted by
        E-value.

    Example:
        >>> queries = ["MAGFLKVVQLLAKYGSKAVQWAWANKGKILDWLNAGQAIDWVVSKIKQILGIK"]
        >>> targets = ([
        ...     "MESILDLQELETSEEESALMAASTVSNNC",
        ...     "MKKAVIVENKGCATCSIGAACLVDGPIPDFEIAGATGLFGLWG",
        ...     "MAGFLKVVQILAKYGSKAVQWAWANKGKILDWINAGQAIDWVVEKIKQILGIK",
        ...     "MTQIKVPTALIASVHGEGQHLFEPMAARCTCTTIISSSSTF",
        ... ])
        >>> for hit in pyswrd.search(queries, targets):
        ...     cigar = hit.result.cigar()
        ...     print(f"target={hit.target_index} score={hit.score} evalue={hit.evalue:.1g} cigar={cigar}")
        target=2 score=268 evalue=1e-33 cigar=53M

    """
    cdef Sequences                      query_db
    cdef Sequences                      target_db
    cdef Scorer                         scorer
    cdef pyopal.lib.ScoreMatrix         score_matrix
    cdef HeuristicFilter                hfilter
    cdef EValue                         evalue
    cdef double                         target_evalue
    cdef size_t                         query_index
    cdef size_t                         target_index
    cdef size_t                         query_length
    cdef size_t                         target_length
    cdef vector[uint32_t]               target_indices
    cdef Sequences                      sub_db
    cdef pyopal.lib.ScoreResult         score_result
    cdef pyopal.lib.FullResult          target_result
    cdef vector[pair[uint32_t, double]] target_evalues

    if threads == 0:
        threads = os.cpu_count() or 1
    if ThreadPool is None:
        threads = 1

    query_db  = queries if isinstance(queries, Sequences) else Sequences(queries)
    target_db = targets if isinstance(targets, Sequences) else Sequences(targets)

    scorer  = Scorer(name=scorer_name, gap_open=gap_open, gap_extend=gap_extend)
    score_matrix = scorer.score_matrix

    with query_db.lock.read, target_db.lock.read:
        # run the filter on the DB
        hfilter = HeuristicFilter(query_db, kmer_length=kmer_length, max_candidates=max_candidates, score_threshold=score_threshold, scorer=scorer, threads=threads)
        filter_result = hfilter.score(target_db).finish()
        evalue = EValue(filter_result.database_length, scorer)
        # align the DB
        with (nullcontext(None) if threads == 1 else ThreadPool(threads)) as pool:
            for query_index, query in enumerate(query_db):
                # get length of query
                query_length = len(query)
                # extract candidates and align them in scoring mode only
                target_indices = filter_result._indices[query_index]
                sub_db = target_db.extract(target_indices)
                score_results = align(query, sub_db, algorithm=algorithm, mode="score", score_matrix=score_matrix, gap_open=gap_open, gap_extend=gap_extend, threads=threads, pool=pool, ordered=True)
                # extract indices with E-value under threshold
                target_evalues.clear()
                for score_result, target_index in zip(score_results, target_indices):
                    target_length = target_db._lengths[target_index]
                    target_evalue = evalue.calculate(score_result.score, query_length, target_length)
                    if target_evalue <= max_evalue:
                        target_evalues.emplace_back(target_index, target_evalue)
                # skip second stage if no E-value passed the threshold
                if target_evalues.empty():
                    continue
                # get only `max_alignments` alignments per query, smallest e-values first
                with nogil:
                    stable_sort_by_second(target_evalues)
                    if target_evalues.size() > max_alignments:
                        target_evalues.resize(max_alignments)
                    target_indices.clear()
                    for target_pair in target_evalues:
                        target_indices.push_back(target_pair.first)
                # align selected sequences
                sub_db = target_db.extract(target_indices)
                ali_results = align(query, sub_db, algorithm=algorithm, mode="full", score_matrix=score_matrix, gap_open=gap_open, gap_extend=gap_extend, threads=threads, pool=pool, ordered=True)
                # return hits for aligned sequences
                for (target_index, target_evalue), target_result in zip(target_evalues, ali_results):
                    yield Hit(
                        query_index=query_index,
                        target_index=target_index,
                        evalue=target_evalue,
                        result=target_result
                    )