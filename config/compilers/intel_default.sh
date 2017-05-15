# Initialize some variables.
omp=
mpi=
mkl=

mpi_libs=
misc_libs=

AR=ar

echo "Intel Fortran Compiler"

MODDIRPREFIX="-module "

# Add some additinal flags to the object directory
if [[ -n $USE_MIC ]]; then
    DPO=${DPO_BASE}/${DPO}_INTEL_MIC/;
else
    DPO=${DPO_BASE}/${DPO}_INTEL/;
fi
if test ! -d $DPO; then  mkdir $DPO; fi

# Set OpenMP flags.
if test $USE_SMP = 1; then omp="-openmp"; fi

# Set compiler commands for DMP or serial.
if test $USE_DMP = 1; then
  FORTRAN_CMD=mpiifort
  LINK_CMD=mpiifort
else
  FORTRAN_CMD=ifort
  LINK_CMD=ifort
fi

# --> Verify compiler is in $PATH <-- #

SET_MPI_INCLUDE
if test $USE_DMP = 1; then
  mpi="-I$MPI_INCLUDE_PATH -DMPI"
  mpi_libs=
fi

# Set generic library information:
ode="${DPO}odepack.a"
blas="${DPO}blas90.a"
dgtsv="${DPO}dgtsv90.a"

mkl_libs="${blas} ${dgtsv}"
LIB_FLAGS="${ode} ${mkl_libs} ${misc_libs} ${mpi_libs} "

# Debug flags for Intel Fortran
dbg=
if test "${USE_DEBUG}" = "1"; then dbg="-g"; fi

# Common compile flags.
common="-c -I. -convert big_endian -assume byterecl -cpp"

if [[ -n $USE_MKL ]]; then
    mkl="-mkl"
    common=${common}" -mkl"
else
    LIB_FLAGS="${LIB_FLAGS} ${blas} ${dgtsv}"
fi

# To display the flags mfix.exe is compiled with, run:
# >
# > readelf -p .debug_str mfix.exe | grep switches

if [[ -n $USE_MIC ]]; then
    common=${common}" -mmic"
else
    common=${common}" -xHost"
fi

if test "${USE_CODECOV}" = "1"; then common=${common}" -prof-gen=srcpos"; fi

case $OPT in
  0)echo "Setting flags for debugging."
    dbg="-traceback -check all -O0"
    if [[ -z $USE_MIC ]]; then dbg="$dbg -fpe:0" ; fi
    FORT_FLAGS="${omp} ${mpi} ${common} -FR ${dbg} -g"
    FORT_FLAGS3="${common} -O0 -g"
    LINK_FLAGS="${omp} ${mkl} -g";;

  1)echo "Setting flags for low optimization."
    FORT_FLAGS="${omp} ${mpi} ${common} -FR -O1 ${dbg}"
    FORT_FLAGS3="${common} -O1 ${dbg}"
    LINK_FLAGS="${omp} ${mkl} ${dbg}";;

  2)echo "Setting flags for medium optimization."
    FORT_FLAGS="${omp} ${mpi} ${common} -FR -O2 ${dbg}"
    FORT_FLAGS3="${common} -O2 ${dbg}"
    LINK_FLAGS="${omp} ${mkl} ${dbg}";;

  3)echo "Setting flags for high optimization."
    AR=xiar
    FORT_FLAGS="${omp} ${mpi} -FR -fast -no-ipo ${common} ${dbg}"
    FORT_FLAGS3="${common} -fast ${dbg}"
    LINK_FLAGS="${omp} ${mkl} ${dbg}";;
  *)echo "Unsupported optimization level."
    exit;;
esac

if [[ -n $USE_MIC ]]; then LINK_FLAGS=${LINK_FLAGS}" -mmic"; fi
