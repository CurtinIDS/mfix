# Initialize some variables.
omp=
mpi=
mkl=

mkl_libs=
mpi_libs=
misc_libs=

echo "HOPPER@NERSC :: Intel Fortran Compiler"

MODDIRPREFIX="-module "

DPO=${DPO_BASE}/${DPO}_NERSC_HOPPER_INTEL/;
if test ! -d $DPO; then  mkdir $DPO; fi

# Set OpenMP flags.
if test $USE_SMP = 1; then omp="-openmp"; fi


# Set compiler commands for DMP or serial.
if test $USE_DMP = 1; then
  FORTRAN_CMD=ftn
  LINK_CMD=ftn
else
  FORTRAN_CMD=ftn
  LINK_CMD=ftn
fi


# --> Verify compiler is in $PATH <-- #


#SET_MPI_INCLUDE
#if test $USE_DMP = 1; then
#  mpi="-I$MPI_INCLUDE_PATH"
#  mpi_libs=
#fi

# Set generic library information:
ode="${DPO}odepack.a"
blas="${DPO}blas90.a"
dgtsv="${DPO}dgtsv90.a"


# The Intel MKL is disabled .
mkl_libs="${blas} ${dgtsv}"


LIB_FLAGS="${ode} ${mkl_libs} ${mpi_libs} ${misc_libs}"

# Debug flags for Intel Fortran
dbg=
if test "${USE_DEBUG}" = "1"; then dbg="-g"; fi

# Common compile flags.
common="-c -I. -convert big_endian -assume byterecl -cpp"

case $OPT in
  0|1|2|3)echo "Setting compiler flags."

    FORT_FLAGS="${omp} ${mpi} ${mkl} ${common} -FR -O3${dbg}"
    FORT_FLAGS3="${common} ${mkl} -O3 ${dbg}"
    LINK_FLAGS="${omp} ${dbg}";;

  *)echo "Unsupported optimization level."
    exit;;
esac
