#if USE_OPENMP
    #include <omp.h>
#endif

// OpenMP loop macros
#if USE_OPENMP
    #define OMP_STATIC_LOOP     _Pragma("omp parallel for schedule(static)")
    #define OMP_DYNAMIC_LOOP    _Pragma("omp parallel for schedule(dynamic, 64)")
#else
    #define OMP_STATIC_LOOP
    #define OMP_DYNAMIC_LOOP
#endif