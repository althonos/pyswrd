from libc.stdint cimport int32_t, uint32_t, uint64_t
from libcpp.memory cimport shared_ptr, unique_ptr

from sword.score_matrix cimport ScoreMatrix

cdef extern from "evalue.hpp" nogil:

    cppclass EValue:
        double calculate(int32_t score, uint32_t query_length, uint32_t target_length)

    unique_ptr[EValue] createEValue(uint64_t database_cells, shared_ptr[ScoreMatrix])