#ifndef AMREX_OPENMP_H_
#define AMREX_OPENMP_H_
#include <AMReX_Config.H>

#ifdef AMREX_USE_OMP
#include <AMReX_Extension.H>
#include <omp.h>
#include <array>

namespace amrex::OpenMP {

    inline int get_num_threads () { return omp_get_num_threads(); }
    inline int get_max_threads () { return omp_get_max_threads(); }
    inline int get_thread_num  () { return omp_get_thread_num();  }
    inline int in_parallel     () { return omp_in_parallel();     }
    inline void set_num_threads (int num) { omp_set_num_threads(num); }

    void Initialize ();
    void Finalize ();

    static constexpr int nlocks = 128;
    extern AMREX_EXPORT std::array<omp_lock_t,nlocks> omp_locks;
}

#else // AMREX_USE_OMP

namespace amrex::OpenMP {

    constexpr int get_num_threads () { return 1; }
    constexpr int get_max_threads () { return 1; }
    constexpr int get_thread_num  () { return 0; }
    constexpr int in_parallel     () { return false; }
    constexpr void set_num_threads (int) { /* nothing */ }
}

#endif // AMREX_USE_OMP

namespace amrex {
    /** ... */
    int
    numUniquePhysicalCores();
}

#endif
