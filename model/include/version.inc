#if ( defined(__GFORTRAN__) &&  ( __GNUC__ > 4 || ( __GNUC__ == 4 && __GNUC_MINOR__ > 5 ) ) ) \
  ||  ( defined(__INTEL_COMPILER) && (__INTEL_COMPILER > 1400) )
#define ERROR_STOP ERROR STOP
#else
#define ERROR_STOP STOP
#endif
