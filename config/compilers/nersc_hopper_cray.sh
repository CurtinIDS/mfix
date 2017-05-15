# Initialize some variables.
omp=
mpi=
mkl=

mkl_libs=
mpi_libs=
misc_libs=

echo "HOPPER@NERSC :: Cray Fortran Compiler"

MODDIRPREFIX="-p "

DPO=${DPO_BASE}/${DPO}_NERSC_HOPPER_CRAY/;
if test ! -d $DPO; then  mkdir $DPO; fi

# Set OpenMP flags.
if test $USE_SMP = 0; then omp="-h noomp"; fi


# Set compiler commands for DMP or serial.
FORTRAN_CMD=ftn
LINK_CMD=ftn

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

# Debug flags for Cray Fortran
dbg=
if test "${USE_DEBUG}" = "1"; then dbg="-rm"; fi

# Common compile flags.
common="-c -h byteswapio -cpp"

case $OPT in
  0|1|2|3)echo "Setting compiler flags."

    FORT_FLAGS="${omp} ${mpi} ${mkl} ${common} -f free ${dbg}"
    FORT_FLAGS3="${common} ${mkl} ${dbg}"
    LINK_FLAGS="${omp} ${dbg}";;

  *)echo "Unsupported optimization level."
    exit;;
esac
