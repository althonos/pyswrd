from libcpp cimport bool

cdef extern from "<mutex>" namespace "std" nogil:
    cdef cppclass mutex:
        mutex() noexcept

        void lock()
        void unlock()
        bool try_lock()

    cdef cppclass unique_lock[Mutex]:
        unique_lock() noexcept
        unique_lock(Mutex& m)

        void lock()
        bool try_lock()
        void unlock()