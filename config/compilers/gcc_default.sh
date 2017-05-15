incs=
mpi=

# Functions for MPI
. $MFIX_CONFIG/mpi_fun.sh

echo "Gfortran on 64 bit machine"

MODDIRPREFIX="-J"

# Add some additinal flags to the object directory
DPO=${DPO_BASE}/${DPO}_GNU/
if test ! -d ${DPO}; then  mkdir ${DPO}; fi

# Add DPO to the includes list.
incs=${incs}" -I${DPO}"

# Set OpenMP flags.
if test $USE_SMP = 1; then
  omp="-fopenmp"
else
  omp=""
fi


# Set MPI flags.
if test $USE_DMP = 1; then
  FORTRAN_CMD=mpif90
  LINK_CMD=mpif90
  mpi="-DMPI"
else
  FORTRAN_CMD=gfortran
  LINK_CMD=gfortran
fi


# Set the MPI include path and add to incs.
SET_MPI_INCLUDE
if test ${USE_DMP} = 1; then
  incs=${incs}" -I${MPI_INCLUDE_PATH}"
fi


# Set library paths:
ode="${DPO}odepack.a"
blas="${DPO}blas90.a"
dgtsv="${DPO}dgtsv90.a"
mpi_libs=
misc_libs=

LIB_FLAGS="${ode} ${blas} ${dgtsv} ${misc_libs} ${mpi_libs} "

# --> Verify compiler is in $PATH <-- #

# Debug flags for GNU Fortran
dbg=
if test ${USE_DEBUG} = 1; then dbg="-g"; fi

# Base flags for GNU Fortran
common="-c -I. ${mpi} ${incs} -fconvert='big-endian' -cpp"
# To display the flags mfix.exe is compiled with, run:
# >
# > readelf -p .GCC.command.line mfix.exe

case $OPT in
  0)echo "Setting flags for debugging."
    dbg="-fbounds-check -fbacktrace -ffpe-trap=invalid,zero,overflow -g"
    FORT_FLAGS="${omp} ${common} -ffree-form -ffree-line-length-0 -O0 ${dbg}"
    FORT_FLAGS3="${omp} ${common} -O0 ${dbg}"
    LINK_FLAGS="${omp} -g";;

  1)echo "Setting flags for low optimization."
    FORT_FLAGS="${omp} ${common} -ffree-form -ffree-line-length-0 -O1 ${dbg}"
    FORT_FLAGS3="${omp} ${common} -O1 ${dbg}"
    LINK_FLAGS="${omp} ${dbg}";;

  2)echo "Setting flags for medium optimization."
    FORT_FLAGS="${omp} ${common} -ffree-form -ffree-line-length-0 -O2 ${dbg}"
    FORT_FLAGS3="${omp} ${common} -O1 ${dbg}"
    LINK_FLAGS="${omp} ${dbg}";;

  3)echo "Setting flags for high optimization."
    FORT_FLAGS="${omp} ${common} -ffree-form -ffree-line-length-0 -O3 ${dbg}"
    FORT_FLAGS3="${omp} ${common} -O2 ${dbg}"
    LINK_FLAGS="${omp} ${dbg}";;

  *)echo "Unsupported optimization level."
    exit;;
esac
